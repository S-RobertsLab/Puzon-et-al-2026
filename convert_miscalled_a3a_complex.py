#!/usr/bin/env python3
"""
correct_snv_ins_complex.py  –  v3 (for wide TSV schema)
-------------------------------------------------------
Right-justify an insertion that comes immediately before a G-based SNV
and (if they end up on the same coordinate) swap the inserted base with
the SNV’s new base, producing e.g.  G>GA   and   G>C.

Works with inputs having columns like (not exhaustive):
CHROM START END CONTEXT MUTATION MUT_TYPE SAMPLE COMPLEX_ID
C_PHASED_READS ... C_CLASS  D_PHASED_READS ... D_CLASS

Notes
-----
• Non-complex rows (COMPLEX_ID in {'.', '', '0'}) are passed through unchanged.
• Only fixes 2-row complexes that are exactly one INS + one SNV.
• Only applies the swap when the INS is a single-base insertion.
• Output order matches the input line-for-line, with each complex emitted
  at its first occurrence in the file.

USAGE
    python correct_snv_ins_complex.py  input.tsv  [-o OUT]

If -o/--output is omitted, <infile>.fixed is written.
"""

import argparse, csv
from pathlib import Path
from collections import defaultdict


# ----------------------------- small helpers ----------------------------- #

def find_upper_idx(s: str) -> int:
    """Return index of the single uppercase letter (assumed unique)."""
    ups = [i for i, ch in enumerate(s) if ch.isupper()]
    return ups[0] if len(ups) == 1 else -1


def move_anchor(ctxt: str, idx: int, shift: int) -> tuple[str, int]:
    """Move the single uppercase base `shift` positions to the right."""
    if shift == 0:
        return ctxt, idx
    clist = list(ctxt)
    clist[idx] = clist[idx].lower()
    new_idx = idx + shift
    clist[new_idx] = clist[new_idx].upper()
    return "".join(clist), new_idx


def right_justify_same_base(ctxt: str, idx: int) -> tuple[str, int]:
    """
    Slide the uppercase base rightward through a homopolymer of *itself*.
    Returns (new_context, extra_shift).
    """
    base = ctxt[idx].lower()
    clist = list(ctxt)
    shift = 0
    pos = idx + 1
    while pos < len(clist) and clist[pos].lower() == base:
        clist[idx + shift] = base          # make previous lowercase
        clist[pos]          = base.upper() # move uppercase
        shift += 1
        pos  += 1
    return "".join(clist), shift


def parse_mut(m: str) -> tuple[str, str]:
    return m.split(">")


def is_single_base_insertion(mut: str) -> bool:
    """len(ALT) == len(REF)+1 for 'REF>ALT'."""
    try:
        ref, alt = mut.split(">")
    except ValueError:
        return False
    return len(alt) == len(ref) + 1


# -------------------------- core logic per complex ------------------------ #

def fix_complex(ins, snv):
    """Return (fixed_INS_row, fixed_SNV_row). Mutates shallow copies (dicts)."""
    s_start = int(snv["START"])
    i_start = int(ins["START"])
    i_end   = int(ins["END"])

    # Only handle INS before SNV and SNV ref == G
    s_ref, s_alt = parse_mut(snv["MUTATION"])
    if i_start >= s_start or s_ref.upper() != "G":
        return ins, snv

    # Only proceed if INS is a single-base insertion
    if not is_single_base_insertion(ins["MUTATION"]):
        return ins, snv

    # Step 1 – move insertion anchor directly onto the SNV coord
    anchor_shift = s_start - i_start          # ≥ 1
    idx = find_upper_idx(ins["CONTEXT"])
    if idx == -1:
        return ins, snv                      # malformed context
    new_ctxt, _ = move_anchor(ins["CONTEXT"], idx, anchor_shift)
    ins["CONTEXT"] = new_ctxt
    ins["START"]   = str(i_start + anchor_shift)
    ins["END"]     = str(i_end   + anchor_shift)

    # Step 2 – homopolymer right-justify on the anchor base
    idx = find_upper_idx(ins["CONTEXT"])
    ins["CONTEXT"], extra = right_justify_same_base(ins["CONTEXT"], idx)
    ins["START"] = str(int(ins["START"]) + extra)
    ins["END"]   = str(int(ins["END"])   + extra)

    # Step 3 – they must now overlap (same START)
    if int(ins["START"]) != s_start:
        return ins, snv

    # Step 4 – swap inserted base ⇄ SNV new base
    _, old_ins_alt = parse_mut(ins["MUTATION"])
    ins_base       = old_ins_alt[-1]         # single-base insertion assumed

    new_ins_base   = s_alt                   # what was SNV alt
    new_snv_base   = ins_base                # what was inserted before swap
    anchor_base    = s_ref.upper()           # “G”

    ins["MUTATION"] = f"{anchor_base}>{anchor_base}{new_ins_base}"
    snv["MUTATION"] = f"{anchor_base}>{new_snv_base}"

    return ins, snv


# ------------------------------- main I/O --------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description="Right-justify INS before G-based SNV and fix bases (wide TSV-safe)."
    )
    ap.add_argument("infile", type=Path)
    ap.add_argument("-o", "--output", type=Path,
                    help="output file (default: <infile>.fixed)")
    args = ap.parse_args()
    if args.output is None:
        args.output = args.infile.with_suffix(".fixed")

    # Read all rows; track complexes and preserve exact file order.
    header = None
    master_rows = []  # rows in file order
    complex_keys_in_order = []  # first occurrence order of complexes
    buckets = defaultdict(list)  # (sample, complex_id) -> [rows in file order]
    first_row_obj = {}  # key -> first row object (to match at emission time)

    def is_simple(cid: str) -> bool:
        cid = (cid or "").strip()
        return cid in (".", "", "0")

    with args.infile.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        header = rdr.fieldnames
        if header is None:
            raise ValueError("Input appears to have no header.")
        required = {"SAMPLE", "COMPLEX_ID", "MUT_TYPE", "MUTATION", "START", "END", "CONTEXT"}
        missing = [c for c in required if c not in header]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")

        for r in rdr:
            master_rows.append(r)
            cid = r.get("COMPLEX_ID", "")
            if is_simple(cid):
                continue  # handled later as pass-through
            key = (r["SAMPLE"], cid)
            if key not in first_row_obj:
                first_row_obj[key] = r
                complex_keys_in_order.append(key)
            buckets[key].append(r)

    # Prepare fixed versions for all eligible complexes.
    fixed_pairs = {}  # key -> [row1_fixed, row2_fixed] in original within-pair order
    for key in complex_keys_in_order:
        rows = buckets[key]
        if len(rows) != 2 or {"INS", "SNV"} != {r["MUT_TYPE"] for r in rows}:
            # Not a 2-member INS+SNV complex → keep as-is
            fixed_pairs[key] = rows
            continue

        # Identify order and run fixer
        ins_row = rows[0] if rows[0]["MUT_TYPE"] == "INS" else rows[1]
        snv_row = rows[1] if rows[0]["MUT_TYPE"] == "INS" else rows[0]

        # Work on shallow copies so we don't mutate shared dicts inadvertently
        ins_copy = dict(ins_row)
        snv_copy = dict(snv_row)
        ins_fixed, snv_fixed = fix_complex(ins_copy, snv_copy)

        # Rebuild in original within-complex order
        if rows[0]["MUT_TYPE"] == "INS":
            fixed_pairs[key] = [ins_fixed, snv_fixed]
        else:
            fixed_pairs[key] = [snv_fixed, ins_fixed]

    # Emit rows in the exact input order.
    out_rows = []
    emitted_complex = set()
    for r in master_rows:
        cid = r.get("COMPLEX_ID", "")
        if is_simple(cid):
            out_rows.append(r)  # pass-through
            continue
        key = (r["SAMPLE"], cid)
        if key in emitted_complex:
            continue  # already emitted this pair at its first occurrence
        # Emit the pair at the first occurrence only.
        if r is first_row_obj.get(key):
            out_rows.extend(fixed_pairs.get(key, buckets[key]))
            emitted_complex.add(key)
        # else: this is the second row of the pair encountered prior to its first?
        # That shouldn't happen, but if it does, just skip (it will be emitted at first).

    with args.output.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    print(f"Wrote {len(out_rows)} rows → {args.output}")


if __name__ == "__main__":
    main()
