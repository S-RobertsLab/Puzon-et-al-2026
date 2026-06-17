#!/usr/bin/env python3
"""
minus_five_to_plus_five_spectra.py
----------------------------------
Flanking-base spectra (percent stacked bars) for SNV, DNV, single-base
DEL (C or G→C) and single-base INS (T or A→T), all oriented to the
pyrimidine strand.

Centres
-------
• DEL  – deleted C
• INS  – inserted T treated as offset 0
• DNV  – two bars ('0a','0b') for the two mutated bases

Optional filtering
------------------
  --apobec-only   keep only C>T and C>G substitutions (SNV) **and**
                  CC→TT/TG/GT/GG substitutions (DNV) **and**
                  single-base C deletions in a tC motif (DEL)

Output: one SVG per mutation class saved beside the input .complex file.
"""

from __future__ import annotations
from pathlib import Path
from collections import defaultdict, Counter
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

##############################################################################
# helpers
##############################################################################
# Fixed base colors (from provided SVG)
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

def rc(seq: str) -> str:
    """reverse complement"""
    return seq.translate(RCMAP)[::-1]

def uidx(s: str) -> list[int]:
    """indices of uppercase letters in s"""
    return [i for i, c in enumerate(s) if c.isupper()]

##############################################################################
#  G-deletion → C-side  (trimmed to essentials)
##############################################################################

def shift_g(ctx: str, g_i: int):
    clist = list(ctx)
    sh, pos = 0, g_i + 1
    while pos < len(clist) and clist[pos].lower() == 'g':
        sh += 1
        pos += 1
    if sh:
        clist[g_i] = clist[g_i].lower()
        clist[g_i + sh] = clist[g_i + sh].upper()
    return "".join(clist), sh

def is_tC_del(ctx: str, mut: str) -> bool:
    """
    Return True iff the event is a single-base C deletion in a tC motif (5' T).
    Assumes the mutated bases in ctx are uppercase; the deleted C is the later
    uppercase letter (ups[1]).
    """
    ups = uidx(ctx)
    if len(ups) != 2 or ups[1] == 0:
        return False
    i2 = ups[1]
    if ctx[i2].upper() != 'C' or ctx[i2 - 1].lower() != 't':
        return False
    ref, alt = mut.split('>')
    return len(ref) == 2 and ref[1] == 'C' and alt == ref[0]

def flip_g_del(start: int, end: int, ctx: str):
    """
    Re-orient a G deletion to the C strand representation.
    If successful and it’s a tC deletion, returns (start, end, flipped_ctx, flipped_mut).
    Otherwise returns None.
    """
    ups = uidx(ctx)
    i2 = ups[1]
    if ctx[i2].upper() != 'G':
        return None
    ctx2, sh = shift_g(ctx, i2)
    start += sh
    end += sh
    rc_ctx = rc(ctx2)
    L = len(rc_ctx)
    g_i = uidx(ctx2)[1]
    c_i = L - 1 - g_i
    lst = list(rc_ctx.lower())
    lst[c_i - 1] = lst[c_i - 1].upper()
    lst[c_i] = lst[c_i].upper()
    flipped = "".join(lst)
    ref = "".join(lst[c_i - 1:c_i + 1]).upper()
    mut = f"{ref}>{ref[0]}"
    if is_tC_del(flipped, mut):
        return start, end, flipped, mut
    return None

def deleted_base(mut: str) -> str:
    """
    For a single-base deletion in minimal VCF-like representation (REF->ALT where
    len(REF) = len(ALT)+1), the removed base is the *last* base of REF.
    Example: 'AC>A' -> removed 'C'.
    """
    ref, alt = mut.split('>')
    if len(ref) != len(alt) + 1:
        raise ValueError("Not a single-base deletion representation")
    return ref[-1].upper()

##############################################################################
# core counting
##############################################################################

def accumulate(path: Path, apobec_only: bool):
    """
    Return flank_counts[mutation_type][offset or '0a'/'0b'][base] = count

    APOBEC rules when apobec_only=True:
      • SNV: C>T or C>G and 5' base == T (tC motif), after orienting to C centre
      • DNV: CC -> TT/TG/GT/GG (both positions mutate to T/G), oriented so first is C
      • DEL: single-base C deletion in a tC motif (if G deletion, re-orient and require tC)
      • INS: no APOBEC constraint (unchanged)
    """
    flank_counts: dict[str, dict[str | int, Counter]] = defaultdict(
        lambda: defaultdict(Counter)
    )

    with path.open() as fh:
        hdr = fh.readline().rstrip().split('\t')
        c = {k: hdr.index(k) for k in
             ("CONTEXT", "MUTATION", "MUT_TYPE", "COMPLEX_ID")}

        for line in fh:
            f = line.rstrip().split('\t')
            # count only simple (non-complex) events
            if f[c["COMPLEX_ID"]] != '.':
                continue

            ctx = f[c["CONTEXT"]]
            mut = f[c["MUTATION"]]
            mtype = f[c["MUT_TYPE"]].upper()

            ########################################################
            # orient & filter
            ########################################################
            if mtype == "SNV":                       # ---------- SNV ----------
                ref, alt = mut.split('>')
                # orient to pyrimidine center
                if ref == 'G':
                    ctx, mut = rc(ctx), rc(mut)
                    ref, alt = mut.split('>')

                # must be C at centre after orientation
                if ref != 'C':
                    continue

                if apobec_only:
                    # APOBEC SNV: C>T or C>G in a tC motif (T immediately 5' of C)
                    if alt not in {'T', 'G'}:
                        continue
                    ups = uidx(ctx)
                    if not ups:
                        continue
                    anchor = ups[0]
                    if anchor == 0 or ctx[anchor - 1].upper() != 'T':
                        continue

            elif mtype == "DNV":                     # ---------- DNV ----------
                ref, alt = mut.split('>')
                if len(ref) != 2 or len(alt) != 2:
                    continue  # skip MNVs / malformed

                # orient so first mutated base is C (pyrimidine-side)
                if ref[0] == 'G':
                    ctx, mut = rc(ctx), rc(mut)
                    ref, alt = mut.split('>')

                if ref[0] != 'C':
                    continue

                if apobec_only:
                    # classic APOBEC DNV: CC -> TT/TG/GT/GG (both Cs mutate to T/G)
                    if ref != 'CC' or any(b not in {'T', 'G'} for b in alt):
                        continue

            elif mtype == "DEL":                     # ---------- DEL ----------
                ref, alt = mut.split('>')
                # single-base deletion only
                if len(ref) != len(alt) + 1:
                    continue

                # determine removed base robustly (no set-difference)
                try:
                    rem = deleted_base(mut)
                except ValueError:
                    continue

                if rem == 'G':
                    # flip to C-side; flip_g_del also enforces tC via is_tC_del()
                    flipped = flip_g_del(0, 0, ctx)
                    if not flipped:
                        continue
                    _, _, ctx, mut = flipped
                    # now it is a C deletion in tC if retained

                elif rem == 'C':
                    # already on C side; if APOBEC-only, require tC motif explicitly
                    if apobec_only and not is_tC_del(ctx, mut):
                        continue
                else:
                    # not a single C/G deletion -> not in scope here
                    continue

            elif mtype == "INS":                     # ---------- INS ----------
                ref, alt = mut.split('>')
                # single-base insertion only
                if len(alt) != len(ref) + 1:
                    continue
                ins = alt[-1].upper()
                # orient inserted base to T (A insertions RC to T)
                if ins == 'A':
                    ctx, mut = rc(ctx), rc(mut)
                elif ins != 'T':
                    continue
                # no APOBEC-only constraint for INS per your spec

            else:
                continue  # unhandled class

            ##############################################
            # count with centring rules
            ##############################################
            ups = uidx(ctx)
            ctxU = ctx.upper()

            if mtype == "DEL":              # centre on deleted C
                anchor = ups[-1]
                for off in range(-5, 6):
                    pos = anchor + off
                    if 0 <= pos < len(ctxU):
                        flank_counts[mtype][off][ctxU[pos]] += 1

            elif mtype == "INS":            # pretend inserted T at 0
                anchor = ups[0]
                for off in range(-5, 0):
                    pos = anchor + off
                    if 0 <= pos < len(ctxU):
                        flank_counts[mtype][off][ctxU[pos]] += 1
                flank_counts[mtype][0]['T'] += 1
                for off in range(1, 6):
                    pos = anchor + off - 1
                    if 0 <= pos < len(ctxU):
                        flank_counts[mtype][off][ctxU[pos]] += 1

            elif mtype == "DNV":            # two bars 0a/0b
                a1, a2 = ups[:2]
                for off in range(-5, 6):
                    lbl = '0a' if off == 0 else off
                    pos = a1 + off
                    if 0 <= pos < len(ctxU):
                        flank_counts[mtype][lbl][ctxU[pos]] += 1
                for off in range(-5, 6):
                    lbl = '0b' if off == 0 else off
                    pos = a2 + off
                    if 0 <= pos < len(ctxU):
                        flank_counts[mtype][lbl][ctxU[pos]] += 1

            else:                           # SNV
                anchor = ups[0]
                for off in range(-5, 6):
                    pos = anchor + off
                    if 0 <= pos < len(ctxU):
                        flank_counts[mtype][off][ctxU[pos]] += 1

    return flank_counts

##############################################################################
# plotting
##############################################################################

def plot(flank_counts, outbase: Path):
    sns.set(style="whitegrid")
    for mtype, posdict in flank_counts.items():
        if mtype == "DNV":
            order = [-5, -4, -3, -2, -1, '0a', '0b', 1, 2, 3, 4, 5]
            xt = [*map(str, [-5, -4, -3, -2, -1]), '0', '0',
                  *map(str, [1, 2, 3, 4, 5])]
        else:
            order = list(range(-5, 6))
            xt = [str(i) for i in order]

        df = (pd.DataFrame(posdict)
              .reindex(index=list("ACGT"))
              .reindex(columns=order)
              .fillna(0).astype(int))

        pct = df.divide(df.sum(axis=0), axis=1) * 100
        # Ensure columns are exactly in A,C,G,T order (they are, from reindex above)
        color_order = [BASE_COLORS[b] for b in pct.T.columns]  # ['A','C','G','T']
        ax = pct.T.plot(kind='bar', stacked=True,
                        color=color_order, edgecolor='black')

        ax.figure.set_size_inches(4, 4)
        ax.set_ylabel('Percentage', fontsize=12, weight='bold')
        ax.set_xlabel('Position', fontsize=12, weight='bold')
        ax.set_title(f'{mtype} ±5 bp base composition', fontsize=12,
                     fontweight='bold')
        ax.set_ylim(0, 100)
        ax.set_xticklabels(xt, rotation=0, fontsize=12, weight='bold')
        h, l = ax.get_legend_handles_labels()
        ax.legend(h[::-1], l[::-1], title='Bases',
            title_fontsize='12', fontsize='12')
        sns.despine()
        svg = outbase.with_name(f"{outbase.stem}_{mtype}_flank.svg")
        ax.figure.tight_layout()
        ax.figure.savefig(svg, format='svg', dpi=300)
        plt.close(ax.figure)

##############################################################################
# entry-point
##############################################################################

def main():
    ap = argparse.ArgumentParser(
        description="±5 bp flanking-base spectra around mutations")
    ap.add_argument('complex_file', help='.complex file')
    ap.add_argument('--apobec-only', action='store_true',
                    help='keep only C>T & C>G SNVs (tC), CC→TT/TG/GT/GG DNVs, and single-base C deletions in tC')
    args = ap.parse_args()

    p = Path(args.complex_file)
    counts = accumulate(p, args.apobec_only)
    plot(counts, p)
    print(f"Plots saved in {p.parent}")

if __name__ == '__main__':
    main()
