#!/usr/bin/env python3
"""
homopolymer_histogram.py
------------------------

Histogram of homopolymer lengths for

  • SNVs :  C>T , C>G , G>A , G>C
  • 1-bp deletions: deletion of *C* or *G*

with an optional --tc-only flag that keeps only
  • TC contexts when the target base is C
  • GA contexts when the target base is G   (reverse-complement rule)

Input columns (tab-separated) must be exactly:

CHROM START END CONTEXT MUTATION MUT_TYPE SAMPLE COMPLEX_ID

Only rows whose COMPLEX_ID == '.' are analysed.
"""

from pathlib import Path
from collections import Counter, defaultdict
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Try to import Savitzky–Golay smoothing; fall back gracefully if unavailable
try:
    from scipy.signal import savgol_filter
    _HAVE_SG = True
except Exception:
    _HAVE_SG = False


# --------------------------------------------------------------------------- #
#                               core utilities                                #
# --------------------------------------------------------------------------- #

def longest_run(text: str, base: str) -> int:
    """Length of the longest uninterrupted run of *base* (already lower-case)."""
    best = cur = 0
    for ch in text:
        if ch == base:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best

def _right_align_idx(context: str, idx: int, base: str) -> int:
    """Slide idx to the rightmost position of a same-base run (case-insensitive)."""
    j = idx
    b = base.lower()
    while (j + 1) < len(context) and context[j + 1].lower() == b:
        j += 1
    return j

def deleted_base_and_pos(before: str, after: str):
    """
    For a single-base deletion MUTATION string 'before>after'
    return (deleted_base_char, index_in_before).
    """
    for i in range(len(before)):
        if i >= len(after) or before[i] != after[i]:
            return before[i], i
    # Fallback: deleted at the end
    return before[-1], len(before) - 1


# --------------------------------------------------------------------------- #
#                       streaming counter  (one file pass)                    #
# --------------------------------------------------------------------------- #

def count_variants(
    file_path: str | Path,
    max_len: int = 20,
    group_by_sample: bool = False,
    tc_only: bool = False,
):
    """
    Returns
    -------
    dict
        { 'SNV': { group: Counter }, 'DEL': { group: Counter } }
    Where each Counter maps homopolymer length (1..max_len) to raw counts.
    Lengths >= max_len are capped into the max_len bin.
    """
    snv_ct: dict[str, Counter] = defaultdict(Counter)
    del_ct: dict[str, Counter] = defaultdict(Counter)

    with Path(file_path).open("rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        needed = ("CONTEXT", "MUTATION", "MUT_TYPE", "SAMPLE", "COMPLEX_ID")
        try:
            col = {name: header.index(name) for name in needed}
        except ValueError as e:
            raise ValueError(f"Input missing required columns {needed}. Got header: {header}") from e

        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[col["COMPLEX_ID"]] != ".":          # skip complex events
                continue

            mut_type = parts[col["MUT_TYPE"]]
            context  = parts[col["CONTEXT"]]
            sample   = parts[col["SAMPLE"]] if group_by_sample else "ALL"

            # -------- indices of the reference allele within CONTEXT --------
            upper_idx = [i for i, c in enumerate(context) if c.isupper()]
            if not upper_idx:
                continue  # malformed line

            if mut_type == "SNV":
                # Contract: exactly ONE uppercase ref base in CONTEXT for SNVs
                if len(upper_idx) != 1:
                    continue

                ref_i   = upper_idx[0]
                ref_base = context[ref_i]          # 'C' or 'G' (we will filter below)
                alt_base = parts[col["MUTATION"]].split(">")[1]

                # keep only C>T/C>G and G>A/G>C
                if (ref_base, alt_base) not in (("C", "T"), ("C", "G"),
                                                ("G", "A"), ("G", "C")):
                    continue

                # TC / GA context filter
                if tc_only:
                    if ref_base == "C":
                        if ref_i == 0 or context[ref_i-1].lower() != "t":
                            continue
                    else:  # ref_base == "G"
                        if ref_i + 1 >= len(context) or context[ref_i+1].lower() != "a":
                            continue

                # APOBEC-style run base:
                #  - If target is C, measure runs of 't'
                #  - If target is G, measure runs of 'a'
                run_base = "t" if ref_base == "C" else "a"

                # Build a "mutated" context for run measurement:
                #  1) make a list so we can overwrite the ref position
                #  2) set the ref position to the run_base (lowercase)
                #  3) lowercase ALL occurrences of run_base (T/A) so contiguous runs
                #     are uniformly lowercase and uninterrupted by capitalization
                mut_ctx_list = list(context)
                mut_ctx_list[ref_i] = run_base
                mut_ctx = "".join(
                    (ch.lower() if ch.lower() == run_base else ch)
                    for ch in mut_ctx_list
                )

                run = min(longest_run(mut_ctx, run_base), max_len)
                snv_ct[sample][run] += 1
                continue

            # -------------------- 1-bp deletion processing -------------------
            if mut_type == "DEL":
                before, after = parts[col["MUTATION"]].split(">")
                if len(before) - len(after) != 1:          # single-bp only
                    continue

                del_base, del_pos_in_before = deleted_base_and_pos(before, after)
                # map that position to CONTEXT index
                if len(upper_idx) != len(before):
                    continue                               # bad formatting
                del_ctx_idx = upper_idx[del_pos_in_before]

                # keep only C or G deletions
                if del_base not in "CG":
                    continue

                # --- Right-align within G-runs for GA checks and for run measurement ---
                aligned_idx = del_ctx_idx
                if del_base == "G":
                    aligned_idx = _right_align_idx(context, del_ctx_idx, "g")

                # TC / GA filter for deletions (after alignment for G)
                if tc_only:
                    if del_base == "C":
                        # require upstream T at the left-aligned C
                        if aligned_idx == 0 or context[aligned_idx - 1].lower() != "t":
                            continue
                    else:
                        # del_base == "G": require A immediately after the rightmost G in the run
                        if (aligned_idx + 1) >= len(context) or context[aligned_idx + 1].lower() != "a":
                            continue

                # choose the run base (t for C-del, a for G-del)
                run_base = "t" if del_base == "C" else "a"

                # build mutated context using the aligned deletion site:
                #   (1) remove the deleted char at aligned_idx
                mut_ctx = context[:aligned_idx] + context[aligned_idx + 1:]
                #   (2) lowercase all run_base letters so runs are continuous
                mut_ctx = mut_ctx.replace(run_base.upper(), run_base)

                run = min(longest_run(mut_ctx, run_base), max_len)
                del_ct[sample][run] += 1

    return {"SNV": snv_ct, "DEL": del_ct}


# --------------------------------------------------------------------------- #
#                                   plotting                                  #
# --------------------------------------------------------------------------- #

def perc(counter: Counter, max_len: int):
    """Return percentage array (len = max_len) summing to 100 over nonzero bins; last bin is ≥ max_len."""
    total = sum(counter.values())
    if total == 0:
        return np.zeros(max_len, dtype=float)
    return np.array([counter.get(i, 0) / total * 100.0 for i in range(1, max_len + 1)])


def _aggregate_counter(group_to_counter: dict[str, Counter]) -> Counter:
    agg = Counter()
    for c in group_to_counter.values():
        agg.update(c)
    return agg


def plot_two_panels(data, file_path, max_len=20):
    """
    Show aggregate % bars (light gray) and overlay per-group % lines.
    This avoids the "non-stacked stacked bars" problem and makes groups comparable.
    """
    types = ("SNV", "DEL")
    x = np.arange(1, max_len + 1)
    fig, axes = plt.subplots(2, 1, figsize=(10, 7.5), sharex=True, sharey=True)

    palette = sns.color_palette("tab10")

    for ax, typ in zip(axes, types):
        counters: dict[str, Counter] = data[typ]
        groups = list(sorted(counters.keys()))

        # ───────── aggregate distribution (bars, % of total) ─────────
        agg = _aggregate_counter(counters)
        y_tot = perc(agg, max_len)  # percentage per bin
        ax.bar(x, y_tot, width=0.85, color="#D9D9D9", edgecolor="#B0B0B0", label="Aggregate")

        # ───────── per-group lines (also % of each group's total) ─────────
        for gi, grp in enumerate(groups):
            y = perc(counters[grp], max_len)
            if y.sum() == 0:
                continue
            ax.plot(x, y, marker="o", lw=1.75, ms=4, label=grp, color=palette[gi % len(palette)])

        # ───────── optional smoothed aggregate curve ─────────
        n_total = sum(agg.values())
        if n_total and _HAVE_SG:
            # pad with a zero at 0 to smooth left edge
            x_curve = np.arange(0, max_len + 1)
            y_curve = np.concatenate([[0.0], y_tot])
            x_fine = np.linspace(0, max_len, max_len * 10 + 1)
            y_interp = np.interp(x_fine, x_curve, y_curve)
            win = 31 if len(x_fine) >= 31 else (len(x_fine) | 1)
            y_smooth = savgol_filter(y_interp, window_length=win, polyorder=3)
            ax.plot(x_fine, y_smooth, color="grey", lw=1.5, alpha=0.8, label="Aggregate (smoothed)")

        # ───────── mean line (capped at max_len) ─────────
        if n_total:
            μ = sum(k * agg[k] for k in range(1, max_len + 1)) / n_total
            ax.axvline(μ, color="red", ls="--", lw=1.8)
            ymax = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1.0
            ax.text(μ + 0.3, ymax * 0.95, f"μ (capped@{max_len}) = {μ:.2f}",
                    color="red", ha="left", va="top", fontsize=9)

        # ───────── cosmetics ─────────
        ax.set_ylabel("Percent of events (%)")
        ax.set_title(f"{typ} (n = {n_total})", fontsize=12)
        ax.grid(axis="y", ls="--", alpha=.5)

        # legend
        if len(groups) > 0:
            ax.legend(frameon=False, fontsize=8, ncol=2)

        # x-axis limits: start at 1, end at max_len
        ax.set_xlim(0.5, max_len + 0.5)

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels([str(i) for i in x])
    axes[-1].set_xlabel("Homopolymer length (bp)")
    sns.despine()
    plt.tight_layout()
    out_svg = file_path.with_name(f"homopolymer_perfection_{file_path.stem}.svg")
    plt.savefig(out_svg, format="svg", dpi=300)
    # plt.show()  # uncomment for interactive preview


# --------------------------------------------------------------------------- #
#                                    CLI                                      #
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(
        description="Homopolymer-length distributions (as %) for APOBEC-related variants"
    )
    ap.add_argument("input", help="variant file (*.complex / TSV)")
    ap.add_argument("--by-sample", action="store_true",
                    help="separate curves for every SAMPLE value")
    ap.add_argument("--max-len", type=int, default=20,
                    help="bucket lengths ≥ this value together")
    ap.add_argument("--tc-only", action="store_true",
                    help="keep only TC (or GA) contexts")
    args = ap.parse_args()

    counts = count_variants(
        args.input,
        max_len=args.max_len,
        group_by_sample=args.by_sample,
        tc_only=args.tc_only,
    )

    input_path = Path(args.input)

    if not any(counts["SNV"].values()):
        print("⚠  No qualifying SNVs found.")
    if not any(counts["DEL"].values()):
        print("⚠  No qualifying deletions found.")

    plot_two_panels(counts, file_path=input_path, max_len=args.max_len)
