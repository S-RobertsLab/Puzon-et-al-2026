#!/usr/bin/env python3
"""
count_complex.py  ·  event counter for *.complex files

Usage:
    python count_complex.py <file.complex>

Creates <file>.counts (TSV) in the same directory.

Counts (per SAMPLE row):
─────────────────────────────────────────────────────────────
• Basic events (only if COMPLEX_ID == '.'):
      SNV, DNV, MNV, INS, DEL
• SNV sub-classes (single-base, non-complex):
      SNV_C_or_G      – ref base C/G
      SNV_A_or_T      – ref base A/T
      SNV_C>T         – C→T or G→A (rev. complement)
      SNV_C>G         – C→G or G→C (rev. complement)
• Single-base deletions (non-complex):
      DEL_1_C_or_G    – deleted base C/G
      DEL_1_A_or_T    – deleted base A/T
• Complex events (COMPLEX_ID ≠ '.'):
      complex_<TYPE1[_TYPE2…]> – one count per unique COMPLEX_ID,
      where TYPE* is any of SNV, DNV, MNV, INS, DEL.
"""

from __future__ import annotations
import argparse
import csv
from collections import defaultdict, Counter
from pathlib import Path
import re
from typing import Tuple
from typing import Dict

def right_justify_indel(
    context: str,
    mutation: str,
    mut_type: str
) -> Tuple[str, int, int, int]:
    """
    Right-justify an insertion or deletion that sits in a perfect repeat and
    give the updated context string.

    Parameters
    ----------
    context   : str   Sequence window with current anchor in UPPER-case.
    mutation  : str   Canonical REF>ALT string (VCF style).
    mut_type  : str   'INS' or 'DEL'.

    Returns
    -------
    new_context : str   Context with the *new* anchor / deleted segment in UPPER-case.
    anchor_idx  : int   0-based index of the new anchor base.
    n_repeats   : int   Number of repeat copies to the right of the original anchor.
    shift_bases : int   Distance (bp) the variant was shifted rightward.
    """
    context = context.strip()

    # ---------------------------------------------------- #
    # 1. Locate current uppercase block (= variant sighting)
    # ---------------------------------------------------- #
    m = re.search(r'[A-Z]+', context)
    if not m:
        raise ValueError("No uppercase letters in context.")
    upper_start, upper_end = m.span()         # inclusive start, exclusive end
    anchor_idx_orig = upper_start
    upper_seg = context[upper_start:upper_end]

    # ---------------------------------------------------- #
    # 2. Extract the sequence that is being added / removed
    # ---------------------------------------------------- #
    if mut_type == "DEL":                     # anchor + deleted bases
        if len(upper_seg) < 2:
            raise ValueError("Deletion must contain ≥1 deleted base.")
        mutated = upper_seg[1:]
    elif mut_type == "INS":
        try:
            ref, alt = mutation.split(">")
        except ValueError:
            raise ValueError("Mutation must be REF>ALT.")
        if not alt.startswith(ref):
            raise ValueError("ALT must begin with REF for an insertion.")
        mutated = alt[len(ref):]              # inserted string
    else:
        raise ValueError("mut_type must be 'INS' or 'DEL'.")

    k = len(mutated)                          # repeat-unit length

    # ---------------------------------------------------- #
    # 3. Count perfect tandem copies to the right
    # ---------------------------------------------------- #
    cursor = anchor_idx_orig + 1 + k          # skip anchor + first unit
    n_repeats = 1                             # already sitting on 1 copy
    while k and context[cursor:cursor+k].upper() == mutated.upper():
        n_repeats += 1
        cursor   += k

    shift_bases = (n_repeats - 1) * k
    anchor_idx  = anchor_idx_orig + shift_bases

    # ---------------------------------------------------- #
    # 4. Build the updated context string
    #    • start with everything lower-case
    #    • upper-case the *new* anchor (+ deleted bases for DEL)
    # ---------------------------------------------------- #
    chars = list(context.lower())
    chars[anchor_idx] = chars[anchor_idx].upper()      # anchor
    if mut_type == "DEL":
        for i in range(k):
            chars[anchor_idx + 1 + i] = chars[anchor_idx + 1 + i].upper()

    new_context = "".join(chars)
    return new_context, anchor_idx, n_repeats, shift_bases

# ---------------------------------------------------------------------- #
#  Reverse complement that preserves the exact upper / lower case pattern
# ---------------------------------------------------------------------- #

# IUPAC bases we want to support
_UPPER = "ACGTRYKMSWBVDHN"
_LOWER = _UPPER.lower()

# Their complements (A↔T, C↔G, R↔Y, K↔M, S↔S, W↔W, B↔V, D↔H, N↔N)
_COMP  = "TGCAYRMKSWVBHDN"
_comp  = _COMP.lower()

# Build one translation table that handles both cases at once
_RC_TRANS: Dict[int, int] = str.maketrans(_UPPER + _LOWER, _COMP + _comp)


def reverse_complement(seq: str) -> str:
    """
    Return the reverse complement of *seq* while preserving the
    lower-/uppercase status of each nucleotide exactly as read in.

    Examples
    --------
    >>> reverse_complement("aCgT")
    'aCgT'
    >>> reverse_complement("tttCTTaaa")
    'tttAAGaaa'
    >>> reverse_complement("GgAaTtCcNn")
    'NnGgAaTtCc'
    """
    # translate → complement; then slice-reverse
    return seq.translate(_RC_TRANS)[::-1]

BASIC_TYPES = {"SNV", "DNV", "MNV", "INS", "DEL"}

def parse_mutation(mutation: str) -> tuple[str, str]:
    ref, alt = mutation.split(">")
    return ref, alt

def deleted_base(ref: str, alt: str) -> str | None:
    """Return the single base deleted if len(ref) - len(alt) == 1, else None"""
    if len(ref) - len(alt) != 1:
        return None
    # simplest heuristic: alt == ref[:-1] (3' deletion) or ref[1:]
    if ref.startswith(alt):
        return ref[len(alt)]
    if ref.endswith(alt):
        return ref[0]
    # fallback: first differing char
    for r, a in zip(ref, alt):
        if r != a:
            return r
    return ref[-1]  # default

def main(complex_path: Path) -> None:
    data = defaultdict(Counter)      # sample → Counter
    complexes = defaultdict(lambda: defaultdict(set))  # sample→ complex_id → set(MUT_TYPE)

    with complex_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample      = row["SAMPLE"]
            mut_type    = row["MUT_TYPE"]
            mutation    = row["MUTATION"]
            complex_id  = row["COMPLEX_ID"]
            context     = row["CONTEXT"]

            if complex_id == ".":
                # Non-complex — count basic categories
                if mut_type in BASIC_TYPES:
                    data[sample][mut_type] += 1

                if mut_type == "SNV":                          # SNV subclasses
                    ref, alt = parse_mutation(mutation)
                    ref_base = ref.upper()
                    alt_base = alt.upper()
                    if "tC" in context or "Ga" in context:
                        data[sample]["SNV_tC"] += 1
                        if mutation in ["C>T", "G>A"]:
                            data[sample]["SNV_tC>T"] += 1
                        if mutation in ["C>G", "G>C"]:
                            data[sample]["SNV_tC>G"] += 1
                    if ref_base in "CG":
                        data[sample]["SNV_C_or_G"] += 1
                    else:
                        data[sample]["SNV_A_or_T"] += 1

                    if (ref_base == "C" and alt_base == "T") or (ref_base == "G" and alt_base == "A"):
                        data[sample]["SNV_C>T"] += 1
                    if (ref_base == "C" and alt_base == "G") or (ref_base == "G" and alt_base == "C"):
                        data[sample]["SNV_C>G"] += 1

                elif mut_type == "DEL":                        # single-base deletions
                    ref, alt = parse_mutation(mutation)
                    del_base = deleted_base(ref, alt)
                    if del_base and len(del_base) == 1:
                        if del_base.upper() in "CG":
                            if "TC" in context:
                                data[sample]["DEL_1_tC"] += 1
                            if del_base == "G":
                                adjusted_context, *_ = right_justify_indel(context, mutation, mut_type)
                                if "Ga" in adjusted_context:
                                    data[sample]["DEL_1_tC"] += 1
                            data[sample]["DEL_1_C_or_G"] += 1
                        else:
                            data[sample]["DEL_1_A_or_T"] += 1
            else:
                # complex; collect composition per complex_id
                complexes[sample][complex_id].add(mut_type)

    # add complex event counts
    for sample, ids in complexes.items():
        for comp_types in ids.values():
            label = "complex_" + "_".join(sorted(comp_types))
            data[sample][label] += 1

    # build column set
    columns = sorted({k for counts in data.values() for k in counts})
    out_path = complex_path.with_suffix(".counts")

    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["SAMPLE", *columns])
        for sample in sorted(data):
            row = [sample] + [data[sample].get(col, 0) for col in columns]
            writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count events in *.complex files")
    parser.add_argument("complex_file", type=Path, help="input .complex TSV")
    args = parser.parse_args()
    main(args.complex_file)
