#!/usr/bin/env python3
"""
filter_snv_ins_complex_phased.py  –  v4 (supports inputs without phasing columns)

Keeps ONLY:
  1) Two–member complexes with exactly one SNV + one single-base INS (len(ALT)=len(REF)+1)
  2) Single (non-complex) DNV rows

Phasing rule (configurable):
  --phasing auto (default): If C_CLASS/D_CLASS exist, use them:
      • Ignore blanks ('.' or '')
      • Reject if ANY class is 'unphased' or 'ambiguous'
      • Accept if AT LEAST ONE class is 'phased'
    If class columns are absent, treat rows as phased (accept).
  --phasing assume-phased: Treat all rows as phased (accept).
  --phasing assume-unphased: Treat all rows as unphased (reject).
  --phasing ignore: Skip phasing checks entirely (accept).

Order is preserved:
  • Each kept complex (both rows) is emitted once at first appearance of either row.
  • Kept single DNV rows are emitted in place.

Output:
  • If -o/--output omitted -> <infile>.snv_ins.phased+dnv.tsv
"""

import argparse
import csv
from pathlib import Path
from collections import defaultdict
from typing import List, Dict, Tuple

SIMPLE_CIDS = {".", "", "0"}

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Keep SNV+1bp-INS complexes and single DNVs with configurable phasing."
    )
    p.add_argument("infile", type=Path, help="Input TSV with header")
    p.add_argument(
        "-o", "--output", type=Path, default=None,
        help="Output TSV (default: <infile>.snv_ins.phased+dnv.tsv)"
    )
    p.add_argument(
        "--phasing",
        choices=["auto", "assume-phased", "assume-unphased", "ignore"],
        default="auto",
        help=(
            "Phasing behavior when class columns (C_CLASS/D_CLASS) are missing or present. "
            "'auto' uses columns if present, otherwise assumes phased."
        ),
    )
    return p.parse_args()

def is_single_base_insertion(mut: str) -> bool:
    """True iff 'REF>ALT' with len(ALT) == len(REF) + 1."""
    try:
        ref, alt = mut.split(">")
    except ValueError:
        return False
    return len(alt) == len(ref) + 1

def row_classes(row: Dict[str, str]) -> List[str]:
    """Collect normalized class labels from C_CLASS and D_CLASS if present in row."""
    vals: List[str] = []
    for col in ("C_CLASS", "D_CLASS"):
        if col in row:
            v = (row.get(col) or "").strip().lower()
            if v and v != ".":
                vals.append(v)
    return vals

def rows_are_phased(
    rows: List[Dict[str, str]],
    phasing_mode: str,
    has_class_cols: bool
) -> bool:
    """
    Decide if rows pass the phasing filter based on mode and column availability.
    """
    if phasing_mode == "ignore":
        return True
    if phasing_mode == "assume-phased":
        return True
    if phasing_mode == "assume-unphased":
        return False

    # phasing_mode == "auto"
    if not has_class_cols:
        # No class columns available → treat as phased (accept)
        return True

    saw_phased = False
    for r in rows:
        for v in row_classes(r):
            if v in {"unphased", "ambiguous"}:
                return False
            if v == "phased":
                saw_phased = True
    return saw_phased

def detect_has_class_cols(header: List[str]) -> bool:
    cols = set(header or [])
    return ("C_CLASS" in cols) or ("D_CLASS" in cols)

def main() -> None:
    args = parse_args()
    if args.output is None:
        args.output = args.infile.with_suffix(".snv_ins.phased+dnv.tsv")

    # --- pass 1: read, preserve order, bucket complexes
    header: List[str] = []
    master_rows: List[Dict[str, str]] = []
    first_row_for_key: Dict[Tuple[str, str], Dict[str, str]] = {}
    order_keys: List[Tuple[str, str]] = []
    buckets: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)

    with args.infile.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        if not header:
            raise ValueError("Input appears to have no header.")

        # Minimal required columns (phasing columns are optional now)
        required = {"SAMPLE", "COMPLEX_ID", "MUT_TYPE", "MUTATION"}
        missing = [c for c in required if c not in header]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")

        has_class_cols = detect_has_class_cols(header)

        for row in reader:
            master_rows.append(row)
            cid = (row.get("COMPLEX_ID") or "").strip()
            if cid in SIMPLE_CIDS:
                continue
            key = (row["SAMPLE"], cid)
            if key not in first_row_for_key:
                first_row_for_key[key] = row
                order_keys.append(key)
            buckets[key].append(row)

    # --- identify which complexes qualify (exactly one INS + one SNV, phased as configured)
    keep_complex = set()
    for key in order_keys:
        rows = buckets[key]
        if len(rows) != 2:
            continue
        types = [r.get("MUT_TYPE", "") for r in rows]
        if sorted(types) != ["INS", "SNV"]:
            continue
        ins_row = rows[0] if rows[0].get("MUT_TYPE") == "INS" else rows[1]
        if not is_single_base_insertion(ins_row.get("MUTATION", "")):
            continue
        if not rows_are_phased(rows, args.phasing, has_class_cols):
            continue
        keep_complex.add(key)

    # --- pass 2: emit in exact input order
    out_rows: List[Dict[str, str]] = []
    emitted_complex = set()
    kept_complex_rows = 0
    kept_dnv_rows = 0

    # If no class cols and phasing='auto', DNV singletons will pass (treated as phased)
    with_class_cols = detect_has_class_cols(header)

    for r in master_rows:
        cid = (r.get("COMPLEX_ID") or "").strip()
        if cid in SIMPLE_CIDS:
            # singleton — only keep if DNV and passes phasing per config
            if r.get("MUT_TYPE") == "DNV":
                if rows_are_phased([r], args.phasing, with_class_cols):
                    out_rows.append(r)
                    kept_dnv_rows += 1
            continue

        key = (r["SAMPLE"], cid)
        if key in emitted_complex:
            continue
        if r is first_row_for_key.get(key):
            if key in keep_complex:
                pair = buckets[key]
                out_rows.extend(pair)
                kept_complex_rows += len(pair)
            emitted_complex.add(key)

    # --- write output
    with args.output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    print(
        f"Kept {kept_complex_rows} rows from SNV+1bp-INS complexes "
        f"({kept_complex_rows//2} complexes) and {kept_dnv_rows} DNV singletons "
        f"(phasing={args.phasing}) → {args.output}"
    )

if __name__ == "__main__":
    main()
