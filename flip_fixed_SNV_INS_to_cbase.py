#!/usr/bin/env python3
"""
flip_g_ins_snv.py  –  v4 (wide TSV-safe)
----------------------------------------
Reverse-complement two-row complexes of the form

        INS (G>Gx)  +  SNV (G>y)

into

        INS (C>Cy)  +  SNV (C>z)

with proper left alignment of the insertion.

• Non-complex rows (COMPLEX_ID in {'.', '', '0'}) are passed through unchanged.
• Only flips 2-row complexes that are exactly one INS + one SNV.
• Guards: SNV must be G-based; INS must be single-base insertion and G>Gx.

USAGE
    python flip_g_ins_snv.py  <corrected.tsv>  [-o OUT]

If -o/--output is not given, <infile>.flipped is written.
"""

import argparse, csv
from pathlib import Path
from collections import defaultdict

COMPL = str.maketrans("ACGTacgt", "TGCAtgca")


# ------------------------------ basic helpers ------------------------------ #

def rc(seq: str) -> str:
    """Reverse complement → lower-case string."""
    return seq.translate(COMPL)[::-1]

def comp(base: str) -> str:
    return base.translate(COMPL).upper()

def parse_mut(m: str):
    ref, alt = m.split(">")
    return ref, alt

def find_upper(s: str) -> int:
    ups = [i for i, ch in enumerate(s) if ch.isupper()]
    return ups[0] if len(ups) == 1 else -1

def is_single_base_insertion(mut: str) -> bool:
    """len(ALT) == len(REF)+1 for 'REF>ALT'."""
    try:
        ref, alt = mut.split(">")
    except ValueError:
        return False
    return len(alt) == len(ref) + 1


# ---------------- left-alignment helper for single-base insertions ---------- #

def left_align_insertion(ctxt: str, idx: int, ins_base: str):
    """
    Slide the *uppercase* anchor (at `idx`) left through a homopolymer of
    `ins_base` (case-insensitive) until the anchor letter differs.
    Returns (new_context, shift_count, new_idx).
    """
    clist = list(ctxt)
    shift = 0
    while idx > 0 and clist[idx - 1].lower() == ins_base.lower():
        clist[idx]     = clist[idx].lower()   # demote current anchor
        idx           -= 1
        clist[idx]     = clist[idx].upper()   # promote previous base
        shift         += 1

    # If anchor still equals insertion base, nudge one more left.
    if idx > 0 and clist[idx].lower() == ins_base.lower():
        clist[idx]     = clist[idx].lower()
        idx           -= 1
        clist[idx]     = clist[idx].upper()
        shift         += 1

    return "".join(clist), shift, idx


# ---------------------------- flip one complex ----------------------------- #

def flip_complex(ins, snv):
    """
    Return (fixed_INS_row, fixed_SNV_row) on the – strand.
    Coordinates refer to +-strand; only START/END of the insertion may shift.
    Assumes: SNV ref == 'G', INS is single-base and G>Gx.
    """
    # Reverse-complement contexts
    ins_ctx_rc = rc(ins["CONTEXT"])
    snv_ctx_rc = rc(snv["CONTEXT"])

    # Add strand annotation
    ins["STRAND"] = snv["STRAND"] = "-"

    # ---- SNV:  G>y  →  C>comp(y)  ----
    s_ref, s_alt      = parse_mut(snv["MUTATION"])
    snv["MUTATION"]   = f"{comp(s_ref)}>{comp(s_alt)}"

    snv_orig_idx      = find_upper(snv["CONTEXT"])
    snv_rc_idx        = len(snv["CONTEXT"]) - 1 - snv_orig_idx
    snv_ctx_list      = list(snv_ctx_rc)
    snv_ctx_list[snv_rc_idx] = snv_ctx_list[snv_rc_idx].upper()
    snv["CONTEXT"]    = "".join(snv_ctx_list)

    # ---- INS:  G>Gx  →  provisional C>Cy (anchor finalised after left-align) ----
    i_ref, i_alt        = parse_mut(ins["MUTATION"])
    ins_base_old        = i_alt[-1]               # x
    ins_base_new        = comp(ins_base_old)      # y
    ins["MUTATION"]     = f"{comp(i_ref)}>{comp(i_ref)}{ins_base_new}"

    ins_orig_idx        = find_upper(ins["CONTEXT"])
    ins_rc_idx          = len(ins["CONTEXT"]) - 1 - ins_orig_idx
    ins_ctx_list        = list(ins_ctx_rc)
    ins_ctx_list[ins_rc_idx] = ins_ctx_list[ins_rc_idx].upper()
    ins_ctx_rc          = "".join(ins_ctx_list)

    # Left-align insertion on RC strand
    ins_ctx_rc, shft, ins_rc_idx = left_align_insertion(
        ins_ctx_rc, ins_rc_idx, ins_base_new
    )
    if shft:
        ins["START"] = str(int(ins["START"]) - shft)
        ins["END"]   = str(int(ins["END"])   - shft)

    # Finalise INS mutation using the (possibly shifted) anchor base
    anchor_base        = ins_ctx_rc[ins_rc_idx].upper()
    ins["CONTEXT"]     = ins_ctx_rc
    ins["MUTATION"]    = f"{anchor_base}>{anchor_base}{ins_base_new}"

    return ins, snv


# ---------------------------------- main ----------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description="Flip G-based INS+SNV complexes to C-based on – strand (wide TSV-safe)."
    )
    ap.add_argument("infile", type=Path)
    ap.add_argument("-o", "--output", type=Path,
                    help="Output (default: <infile>.flipped)")
    args = ap.parse_args()
    if args.output is None:
        args.output = args.infile.with_suffix(".flipped")

    # Read; preserve exact file order; bucket complexes
    header = None
    master_rows = []                  # rows in file order
    complex_keys_in_order = []        # first occurrence order of complexes
    first_row_obj = {}                # key -> first row object
    buckets = defaultdict(list)       # (sample, complex_id) -> rows

    def is_simple(cid: str) -> bool:
        cid = (cid or "").strip()
        return cid in (".", "", "0")

    with args.infile.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        header = rdr.fieldnames
        if header is None:
            raise ValueError("Input appears to have no header.")
        # ensure STRAND exists in output
        if "STRAND" not in header:
            header = header + ["STRAND"]

        required = {"SAMPLE", "COMPLEX_ID", "MUT_TYPE", "MUTATION", "START", "END", "CONTEXT"}
        missing = [c for c in required if c not in rdr.fieldnames]
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(missing)}")

        for row in rdr:
            master_rows.append(row)
            cid = row.get("COMPLEX_ID", "")
            if is_simple(cid):
                continue
            key = (row["SAMPLE"], cid)
            if key not in first_row_obj:
                first_row_obj[key] = row
                complex_keys_in_order.append(key)
            buckets[key].append(row)

    # Prepare flipped versions for eligible complexes
    fixed_pairs = {}  # key -> [row1_fixed, row2_fixed] in original within-pair order
    for key in complex_keys_in_order:
        rows = buckets[key]
        if len(rows) != 2 or {"INS", "SNV"} != {r["MUT_TYPE"] for r in rows}:
            fixed_pairs[key] = rows
            continue

        # Identify order
        ins_row = rows[0] if rows[0]["MUT_TYPE"] == "INS" else rows[1]
        snv_row = rows[1] if rows[0]["MUT_TYPE"] == "INS" else rows[0]

        # Guards: SNV must be G-based; INS must be single-base and G>Gx
        s_ref, _ = parse_mut(snv_row["MUTATION"])
        i_ref, i_alt = parse_mut(ins_row["MUTATION"])
        if s_ref.upper() != "G" or i_ref.upper() != "G" or not is_single_base_insertion(ins_row["MUTATION"]):
            fixed_pairs[key] = rows  # not our pattern → pass through unchanged
            continue

        # Work on shallow copies; flip
        ins_copy = dict(ins_row)
        snv_copy = dict(snv_row)
        ins_fixed, snv_fixed = flip_complex(ins_copy, snv_copy)

        # Rebuild in original within-complex order
        if rows[0]["MUT_TYPE"] == "INS":
            fixed_pairs[key] = [ins_fixed, snv_fixed]
        else:
            fixed_pairs[key] = [snv_fixed, ins_fixed]

    # Emit rows in exact input order; emit each complex once (at first row)
    out_rows = []
    emitted_complex = set()
    for r in master_rows:
        cid = r.get("COMPLEX_ID", "")
        if is_simple(cid):
            # Ensure STRAND field exists in pass-through rows
            if "STRAND" not in r:
                r = dict(r)
                r["STRAND"] = ""
            out_rows.append(r)
            continue

        key = (r["SAMPLE"], cid)
        if key in emitted_complex:
            continue
        # Emit pair at the first occurrence only
        if r is first_row_obj.get(key):
            pair = fixed_pairs.get(key, buckets[key])
            # ensure STRAND field exists on both rows
            fixed_pair = []
            for pr in pair:
                if "STRAND" not in pr:
                    pr = dict(pr)
                    pr["STRAND"] = ""
                fixed_pair.append(pr)
            out_rows.extend(fixed_pair)
            emitted_complex.add(key)

    with args.output.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    print(f"Wrote {len(out_rows)} rows → {args.output}")


if __name__ == "__main__":
    main()
