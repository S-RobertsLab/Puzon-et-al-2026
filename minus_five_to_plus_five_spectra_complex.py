#!/usr/bin/env python3
"""
dnv_flank_phased_singleton.py
-----------------------------
Flanking-base spectra (percent stacked bars) for **DNV only**.

Keep ONLY rows that satisfy ALL of:
  • MUT_TYPE == 'DNV'
  • Not in a complex: COMPLEX_ID is '.' or '' (anything else is ignored)
  • Phased by D_CLASS ONLY: D_CLASS == 'phased' (case-insensitive)
    (C_* columns are completely ignored)

Orientation:
  • Reverse-complement (context + mutation) so the first mutated base is C.
    I.e., after orientation, ref[0] == 'C'.

APOBEC mode (--apobec-only):
  • Keep only classic APOBEC DNVs: CC -> TT/TG/GT/GG (after orientation).

Output:
  <input_stem>_DNV_flank.svg
"""

from __future__ import annotations
from pathlib import Path
from collections import defaultdict, Counter
import argparse
import csv
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Fixed base colors
BASE_COLORS = {
    'A': '#a8c7d6',
    'C': '#8e24aa',
    'G': '#c5d3a3',
    'T': '#d4b7c4',
}

RCMAP = str.maketrans(
    "ACGTacgtRYMKSWHBVDNrymkswhbvdn",
    "TGCAtgcaYRKMWSHVBDByrkmwshvbdn"
)

SIMPLE_CIDS = {".", ""}

def rc(seq: str) -> str:
    """reverse complement"""
    return seq.translate(RCMAP)[::-1]

def uidx(s: str) -> list[int]:
    """indices of uppercase letters in s"""
    return [i for i, c in enumerate(s) if c.isupper()]

def d_class_is_phased(row: dict) -> bool:
    """Only D_CLASS matters; return True iff it's exactly 'phased' (case-insensitive)."""
    v = (row.get("D_CLASS") or "").strip().lower()
    return v == "phased"

def is_singleton_complex_id(row: dict) -> bool:
    """Treat '.' or '' as singleton (non-complex). Anything else is considered a complex."""
    cid = (row.get("COMPLEX_ID") or "").strip()
    return cid in SIMPLE_CIDS

def accumulate(path: Path, apobec_only: bool):
    """
    Return flank_counts['DNV'][offset or '0a'/'0b'][base] = count
    using only phased-singleton DNVs (per rules above).
    """
    flank_counts: dict[str, dict[str | int, Counter]] = defaultdict(
        lambda: defaultdict(Counter)
    )

    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []

        # sanity: required columns
        needed = {"CONTEXT", "MUTATION", "MUT_TYPE", "COMPLEX_ID", "D_CLASS"}
        missing = [c for c in needed if c not in header]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        for r in reader:
            # Skip comment/blank lines gracefully if present
            if not r or all(v == "" or v is None for v in r.values()):
                continue

            if r.get("MUT_TYPE", "").upper() != "DNV":
                continue
            if not is_singleton_complex_id(r):
                continue  # must NOT be in a complex
            if not d_class_is_phased(r):
                continue  # must be phased by D_CLASS only

            ctx = r["CONTEXT"]
            mut = r["MUTATION"]
            if ">" not in mut:
                continue
            ref, alt = mut.split(">")

            # Must be exactly 2→2 bases
            if len(ref) != 2 or len(alt) != 2:
                continue

            # Orient so first mutated base is C
            if ref[0].upper() == "G":
                ctx, mut = rc(ctx), rc(mut)
                ref, alt = mut.split(">")

            if ref[0].upper() != "C":
                continue

            # APOBEC filter if requested
            if apobec_only:
                if ref.upper() != "CC" or any(b.upper() not in {"T", "G"} for b in alt):
                    continue

            # Count ±5 around both mutated positions (two bars: '0a' and '0b')
            ups = uidx(ctx)
            if len(ups) < 2:
                continue
            a1, a2 = ups[:2]
            ctxU = ctx.upper()
            L = len(ctxU)

            for off in range(-5, 6):
                lbl = '0a' if off == 0 else off
                pos = a1 + off
                if 0 <= pos < L:
                    b = ctxU[pos]
                    if b in "ACGT":
                        flank_counts["DNV"][lbl][b] += 1

            for off in range(-5, 6):
                lbl = '0b' if off == 0 else off
                pos = a2 + off
                if 0 <= pos < L:
                    b = ctxU[pos]
                    if b in "ACGT":
                        flank_counts["DNV"][lbl][b] += 1

    return flank_counts

def plot(flank_counts, outbase: Path):
    sns.set(style="whitegrid")

    # Only DNV will be present
    posdict = flank_counts.get("DNV", {})
    order = [-5, -4, -3, -2, -1, '0a', '0b', 1, 2, 3, 4, 5]
    xt = [*map(str, [-5, -4, -3, -2, -1]), '0', '0', *map(str, [1, 2, 3, 4, 5])]

    df = (pd.DataFrame(posdict)
          .reindex(index=list("ACGT"))
          .reindex(columns=order)
          .fillna(0).astype(int))

    # Avoid div-by-zero
    col_sums = df.sum(axis=0).replace(0, pd.NA)
    pct = (df.divide(col_sums, axis=1) * 100).fillna(0)

    color_order = [BASE_COLORS[b] for b in pct.T.columns]  # ['A','C','G','T']
    ax = pct.T.plot(kind='bar', stacked=True, color=color_order, edgecolor='black')

    ax.figure.set_size_inches(4, 4)
    ax.set_ylabel('Percentage', fontsize=12, weight='bold')
    ax.set_xlabel('Position', fontsize=12, weight='bold')
    ax.set_title('DNV ±5 bp base composition (phased singletons)', fontsize=12,
                 fontweight='bold')
    ax.set_ylim(0, 100)
    ax.set_xticklabels(xt, rotation=0, fontsize=12, weight='bold')
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[::-1], l[::-1], title='Bases', title_fontsize='12', fontsize='12')
    sns.despine()

    svg = outbase.with_name(f"{outbase.stem}_DNV_flank.svg")
    ax.figure.tight_layout()
    ax.figure.savefig(svg, format='svg', dpi=300)
    plt.close(ax.figure)

def main():
    ap = argparse.ArgumentParser(
        description="±5 bp flanking-base spectra for phased, non-complex DNVs (wide TSV)."
    )
    ap.add_argument('tsv', help='Wide TSV with CONTEXT/MUTATION/MUT_TYPE/COMPLEX_ID/D_CLASS')
    ap.add_argument('--apobec-only', action='store_true',
                    help='Restrict DNVs to CC→TT/TG/GT/GG after orientation')
    args = ap.parse_args()

    p = Path(args.tsv)
    counts = accumulate(p, args.apobec_only)
    plot(counts, p)
    print(f"DNV plot saved in {p.parent}")

if __name__ == '__main__':
    main()
