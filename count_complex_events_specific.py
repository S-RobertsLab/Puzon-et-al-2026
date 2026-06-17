#!/usr/bin/env python3
"""
count_complex_specific.py  –  counts basic, DNV (COSMIC-78) and
                              user-defined SNV-pair complex events.

Now evaluates SNV-pairs in two biologically relevant views:
  • native (a, b)
  • reverse-complement with reversed order (RC(b), RC(a))

Per-complex de-duplication: a given label is counted at most once across views.
Different labels from different views are both counted.

NEW:
  --annotate            Write an annotated copy of the input with a CLASS_LABELS column.
                        Default OFF.
  --annotate-out PATH   Optional custom path for the annotated TSV (default: <file>.classified.tsv)

Usage:
    python count_complex_specific.py <file.complex> [--annotate [--annotate-out PATH]]
Creates:
  - <file>.specific_counts (TSV)  – per-sample counts
  - <file>.classified.tsv (TSV)   – original rows + CLASS_LABELS (only when --annotate)
"""

from __future__ import annotations
import argparse, csv, re
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, Tuple, List, Set, Optional

###############################################################################
#  Reverse-complement that preserves upper/lower case
###############################################################################
_UPPER = "ACGTRYKMSWBVDHN"
_LOWER = _UPPER.lower()
_COMP  = "TGCAYRMKSWVBHDN"
_comp  = _COMP.lower()
_RC_TRANS: Dict[int, int] = str.maketrans(_UPPER + _LOWER, _COMP + _comp)

def rc(seq: str) -> str:
    "Reverse complement, preserving case."
    return seq.translate(_RC_TRANS)[::-1]

###############################################################################
#  IUPAC degeneracy helper + complements
###############################################################################
IUPAC_CODES = {
    "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"),
    "R": set("AG"), "Y": set("CT"), "W": set("AT"), "S": set("CG"),
    "M": set("AC"), "K": set("GT"), "D": set("AGT"), "H": set("ACT"),
    "B": set("CGT"), "V": set("ACG"), "N": set("ACGT"),
}

IUPAC_COMP = {
    "A":"T","C":"G","G":"C","T":"A",
    "R":"Y","Y":"R","W":"W","S":"S",
    "M":"K","K":"M","D":"H","H":"D",
    "B":"V","V":"B","N":"N"
}

def _comp_code(code: str) -> str:
    code = code.upper()
    if code not in IUPAC_COMP:
        raise ValueError(f"Unknown IUPAC code: {code}")
    return IUPAC_COMP[code]

def rc_code_pair(patt: str) -> str:
    """Reverse-complement a pattern like 'C>D' (IUPAC-aware)."""
    ref, alt = patt.split(">")
    return f"{_comp_code(ref)}>{_comp_code(alt)}"

def rc_snv(snv: str) -> str:
    """Reverse-complement a concrete/IUPAC SNV like 'G>A'."""
    ref, alt = snv.split(">")
    return f"{_comp_code(ref)}>{_comp_code(alt)}"

###############################################################################
#  DNV canonicalisation (COSMIC 78-class – strand-agnostic, pyrimidine-rich)
###############################################################################
def canonicalise_dnv(mutation: str) -> str:
    """Return COSMIC-style canonical dinucleotide change, e.g. CT>AA."""
    ref, alt = mutation.upper().split(">")
    assert len(ref) == len(alt) == 2, f"DNV must be 2 bp: {mutation}"
    pyris = lambda s: s.count("C") + s.count("T")
    ref_f, alt_f = ref, alt
    ref_r, alt_r = rc(ref), rc(alt)
    if pyris(ref_r) > pyris(ref_f) or (pyris(ref_r) == pyris(ref_f) and ref_r < ref_f):
        ref_f, alt_f = ref_r, alt_r
    return f"{ref_f}>{alt_f}"

###############################################################################
#  SNV-pair complex classification helpers
###############################################################################
PairLabel = str

def _match_snv(snv: str, pattern: str) -> bool:
    """Does a mutation (e.g. 'C>T') satisfy a degeneracy pattern (e.g. 'C>D')?"""
    ref, alt = snv.split(">")
    patt_ref, patt_alt = pattern.split(">")
    return (ref.upper() in IUPAC_CODES[patt_ref] and
            alt.upper() in IUPAC_CODES[patt_alt])

def _both(a: str, b: str, p: str) -> bool:
    """Both SNVs match the same pattern p."""
    return _match_snv(a, p) and _match_snv(b, p)

def _unordered(a: str, b: str, p1: str, p2: str) -> bool:
    """Order-insensitive two-pattern match."""
    return (_match_snv(a, p1) and _match_snv(b, p2)) or \
           (_match_snv(a, p2) and _match_snv(b, p1))

def _ordered(a: str, b: str, p1: str, p2: str) -> bool:
    """Order-sensitive two-pattern match (a matches p1 and b matches p2)."""
    return _match_snv(a, p1) and _match_snv(b, p2)

def classify_snv_pair_union(snvs: List[str]) -> Set[PairLabel]:
    """
    Return the **set** of matching labels across exactly two views:
      V0 = (a, b)
      V1 = (RC(b), RC(a))

    De-duplicates per label across views (i.e., a label appears at most once).
    Different labels from different views are both included.
    """
    a, b = snvs
    labels: Set[PairLabel] = set()

    views = [
        ("native",   a,          b),
        ("rc_swap",  rc_snv(b),  rc_snv(a)),
    ]

    for _, x, y in views:
        # Unordered classes
        if _both(x, y, "C>T"):                 labels.add("SNVpair_CT_CT")
        if _both(x, y, "C>D"):                 labels.add("SNVpair_CD_CD")
        if _unordered(x, y, "C>T", "G>A"):     labels.add("SNVpair_CT_GA")
        if _unordered(x, y, "C>R", "G>Y"):     labels.add("SNVpair_CR_GY")
        if _both(x, y, "W>N"):                 labels.add("SNVpair_WN_WN")

        # Ordered classes (direction matters; distinct labels kept)
        # if _ordered(x, y, "C>D", "D>N"):       labels.add("SNVpair_CD_then_DN")
        # if _ordered(x, y, "D>N", "C>D"):       labels.add("SNVpair_DN_then_CD")
        if _ordered(x, y, "A>N", "C>D"):       labels.add("SNVpair_DN_then_CD")
        if _ordered(x, y, "T>N", "C>D"):       labels.add("SNVpair_DN_then_CD")
        if _ordered(x, y, "G>Y", "C>D"):       labels.add("SNVpair_DN_then_CD")

        if _ordered(x, y, "C>D", "A>N"):       labels.add("SNVpair_CD_then_DN")
        if _ordered(x, y, "C>D", "T>N"):       labels.add("SNVpair_CD_then_DN")
        if _ordered(x, y, "C>D", "G>Y"):       labels.add("SNVpair_CD_then_DN")
        
        return labels

def _is_tC_context(row: dict) -> bool:
    """
    True iff the SNV is a cytosine preceded by thymine on the plus
    strand (tC) **or** an adenine preceded by guanine on the minus
    strand (gA). Context string has exactly one uppercase ref base.
    """
    ctx = row.get("CONTEXT")
    if not ctx:
        return False
    m = re.search(r'[ACGT]', ctx)
    if not m:
        return False
    i = m.start()
    if i == 0 or i == len(ctx) - 1:
        return False
    tri = ctx[i-1:i+2]
    return ("tC" in tri) or ("Ga" in tri)

###############################################################################
#  Main counter + optional annotation
###############################################################################
BASIC_TYPES = {"SNV", "DNV"}

def main(complex_file: Path, annotate: bool, annotate_out: Optional[Path]):
    data       = defaultdict(Counter)                  # SAMPLE → Counter
    complexes  = defaultdict(lambda: defaultdict(list))# SAMPLE→ID→row-dicts
    all_rows: List[dict] = []                          # preserve original rows in order
    header_fields: List[str] = []

    # ────────────────────────────────────────────────────────────────
    #  First pass – collect rows & basic (simple) counts
    # ────────────────────────────────────────────────────────────────
    with complex_file.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header_fields = reader.fieldnames or []
        for idx, row in enumerate(reader, start=0):
            # Keep a copy with file order for later annotation output
            row_copy = dict(row)
            row_copy["_idx"] = idx
            all_rows.append(row_copy)

            sample, mtype, mutation, cid = \
                row["SAMPLE"], row["MUT_TYPE"], row["MUTATION"], row["COMPLEX_ID"]

            if cid == ".":                              # simple, non-complex
                if mtype in BASIC_TYPES:
                    data[sample][mtype] += 1
                if mtype == "DNV":                      # 78-class split
                    label = f"DNV_{canonicalise_dnv(mutation)}"
                    data[sample][label] += 1
            else:                                       # part of a complex
                complexes[sample][cid].append(row | {"_idx": idx})  # keep order

    # ────────────────────────────────────────────────────────────────
    #  Second pass – analyse 2-SNV complexes (two views, de-duped labels)
    #  Also build per-complex label sets for optional annotation.
    # ────────────────────────────────────────────────────────────────
    complex_labels: Dict[Tuple[str, str], Set[str]] = {}

    for sample, cid2rows in complexes.items():
        for cid, rows in cid2rows.items():
            snv_rows = [r for r in rows if r["MUT_TYPE"] == "SNV"]
            if len(snv_rows) != 2:
                # Mark non-2-SNV complexes for annotation visibility
                complex_labels[(sample, cid)] = {"NON_2SNV_COMPLEX"}
                continue

            snv_rows.sort(key=lambda r: r["_idx"])      # preserve file order
            a_b = [r["MUTATION"] for r in snv_rows]

            # Primary labels: union across (a,b) and (RC(b),RC(a))
            labels = classify_snv_pair_union(a_b)

            # Secondary (independent) label: tC/tC (gA/gA) on original contexts
            if all(_is_tC_context(r) for r in snv_rows):
                labels = set(labels)  # ensure it's a set
                labels.add("SNVpair_tC_tC")

            # Record labels for counting
            for lab in labels:
                data[sample][lab] += 1

            # Save labels for annotation
            complex_labels[(sample, cid)] = labels

    # ────────────────────────────────────────────────────────────────
    #  Output counts
    # ────────────────────────────────────────────────────────────────
    columns = sorted({k for counts in data.values() for k in counts})
    out_counts = complex_file.with_suffix(".specific_counts")

    with out_counts.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["SAMPLE", *columns])
        for sample in sorted(data):
            writer.writerow([sample] + [data[sample].get(c, 0) for c in columns])

    # ────────────────────────────────────────────────────────────────
    #  Optional annotation output
    # ────────────────────────────────────────────────────────────────
    if annotate:
        out_anno = annotate_out if annotate_out is not None else complex_file.with_suffix(".classified.tsv")
        # Build a lookup from (SAMPLE, COMPLEX_ID) → labels
        # For simple rows, also derive per-row labels where applicable.
        class_col = "CLASS_LABELS"
        final_header = header_fields + [class_col] if class_col not in header_fields else header_fields

        with out_anno.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=final_header)
            writer.writeheader()

            for row in sorted(all_rows, key=lambda r: r["_idx"]):
                sample = row["SAMPLE"]
                cid    = row["COMPLEX_ID"]
                mtype  = row["MUT_TYPE"]
                mutation = row["MUTATION"]

                labels_out: List[str] = []

                if cid != ".":
                    labels = complex_labels.get((sample, cid), set())
                    if labels:
                        labels_out = sorted(labels)
                else:
                    # Simple rows: annotate DNV canonical class; leave SNV blank
                    if mtype == "DNV":
                        labels_out = [f"DNV_{canonicalise_dnv(mutation)}"]

                out_row = dict(row)
                out_row.pop("_idx", None)
                out_row[class_col] = ";".join(labels_out) if labels_out else ""
                writer.writerow(out_row)

###############################################################################
#  CLI
###############################################################################
if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Count events in *.complex files with extra SNV-pair and DNV detail; optional annotated copy."
    )
    p.add_argument("complex_file", type=Path)
    p.add_argument(
        "--annotate",
        action="store_true",
        help="Also write an annotated copy of the input with a CLASS_LABELS column (default: off)."
    )
    p.add_argument(
        "--annotate-out",
        type=Path,
        default=None,
        help="Path for annotated TSV (default: <file>.classified.tsv)."
    )
    args = p.parse_args()
    main(args.complex_file, args.annotate, args.annotate_out)
