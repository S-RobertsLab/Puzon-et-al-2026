#!/usr/bin/env python3
"""
classify_homopolymers.py  •  v2 (wide TSV-safe) •  2025-06-30
--------------------------------------------------------------
Annotate each INS+SNV *two-row* complex with a homopolymer category and
prepend a commented count table. Non-complex or non-INS+SNV rows are
passed through with CATEGORY='other' but are *not* included in the counts.

USAGE
    python classify_homopolymers.py  <in.tsv>  [-o OUT]

If OUT is omitted, the script writes <infile>.classified.tsv next to <in.tsv>.
"""

import argparse, csv
from pathlib import Path
from collections import defaultdict, Counter

# ────────────────────────────────────────────────────────────────────────────
# helpers
# ────────────────────────────────────────────────────────────────────────────

def find_upper(s: str) -> int:
    ups = [i for i, ch in enumerate(s) if ch.isupper()]
    return ups[0] if len(ups) == 1 else -1

def parse_mut(m: str):
    ref, alt = m.split(">")
    return ref, alt

def is_simple_complex_id(cid: str) -> bool:
    cid = (cid or "").strip()
    return cid in (".", "", "0")

# ────────────────────────────────────────────────────────────────────────────
# classification logic (same semantics as your original)
# ────────────────────────────────────────────────────────────────────────────

def classify_complex(ins, snv) -> str:
    """Return one of the seven category strings, or 'other'."""
    if parse_mut(snv["MUTATION"])[0].upper() != "C":
        return "other"

    ins_start, snv_start = int(ins["START"]), int(snv["START"])
    if ins_start == snv_start:
        return "other"

    five_prime = ins_start < snv_start   # insertion 5′ of C?
    orient = "5′" if five_prime else "3′"

    idx_ins = find_upper(ins["CONTEXT"])
    idx_c   = find_upper(snv["CONTEXT"])
    if idx_ins == -1 or idx_c == -1:
        return "other"

    # Map insertion anchor into the SNV context string
    ctx = snv["CONTEXT"]
    offset = snv_start - ins_start if five_prime else ins_start - snv_start
    idx_ins_ctx = idx_c - offset if five_prime else idx_c + offset
    if not (0 <= idx_ins_ctx < len(ctx)):
        return "other"

    ins_base = parse_mut(ins["MUTATION"])[1][-1].lower()

    if five_prime:
        between  = ctx[idx_ins_ctx + 1 : idx_c]
        adjacent = ctx[idx_c - 1].lower() if idx_c - 1 >= 0 else ""
    else:
        between  = ctx[idx_c + 1 : idx_ins_ctx]
        adjacent = ctx[idx_c + 1].lower() if idx_c + 1 < len(ctx) else ""

    if adjacent == "":
        return "other"

    connected = adjacent == ins_base and all(ch.lower() == ins_base for ch in between)

    if ins_base == 't':
        return f"{orient} {'connected' if connected else 'disconnected'} T homopolymer"
    else:
        return f"{orient} disconnected D homopolymer"

# ────────────────────────────────────────────────────────────────────────────
# main
# ────────────────────────────────────────────────────────────────────────────

def main():
    pa = argparse.ArgumentParser(description="Classify INS+SNV homopolymer context (wide TSV-safe).")
    pa.add_argument("infile", type=Path)
    pa.add_argument("-o", "--output", type=Path,
                    help="Output TSV (default: <infile>.classified.tsv)")
    args = pa.parse_args()

    # default → same folder, original name + '.classified.tsv'
    if args.output is None:
        args.output = args.infile.with_name(args.infile.name + ".classified.tsv")

    # Read preserving exact file order; bucket true complexes by (SAMPLE, COMPLEX_ID)
    header = None
    master_rows = []                     # rows in file order
    complex_first_row = {}               # key -> first row object
    complex_keys_in_order = []           # order of first occurrence
    buckets = defaultdict(list)          # key -> rows in file order

    with args.infile.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        if rdr.fieldnames is None:
            raise ValueError("Input appears to have no header.")
        header = list(rdr.fieldnames)
        for row in rdr:
            master_rows.append(row)
            cid = row.get("COMPLEX_ID", "")
            if is_simple_complex_id(cid):
                continue  # singleton/simple rows will be handled pass-through
            key = (row["SAMPLE"], row["COMPLEX_ID"])
            if key not in complex_first_row:
                complex_first_row[key] = row
                complex_keys_in_order.append(key)
            buckets[key].append(row)

    # Ensure CATEGORY exists in output header
    if "CATEGORY" not in header:
        header = header + ["CATEGORY"]

    # Compute category per *two-row INS+SNV* complex; others -> 'other' but not counted
    per_complex_category = {}  # key -> category
    counts = Counter()

    for key in complex_keys_in_order:
        rows = buckets[key]
        if len(rows) == 2 and {"INS", "SNV"} == {r["MUT_TYPE"] for r in rows}:
            ins = rows[0] if rows[0]["MUT_TYPE"] == "INS" else rows[1]
            snv = rows[1] if rows[0]["MUT_TYPE"] == "INS" else rows[0]
            cat = classify_complex(ins, snv)
            per_complex_category[key] = cat
            counts[cat] += 1
        else:
            per_complex_category[key] = "other"  # not counted

    ordered = ["5′ connected T homopolymer",
               "5′ disconnected T homopolymer",
               "3′ connected T homopolymer",
               "3′ disconnected T homopolymer",
               "5′ disconnected D homopolymer",
               "3′ disconnected D homopolymer",
               "other"]

    # Emit rows in exact input order, annotating CATEGORY
    out_rows = []
    emitted_complex = set()

    for r in master_rows:
        cid = r.get("COMPLEX_ID", "")
        if is_simple_complex_id(cid):
            # pass-through singleton rows; label as 'other'
            rr = dict(r)
            rr["CATEGORY"] = rr.get("CATEGORY", "other") or "other"
            out_rows.append(rr)
            continue

        key = (r["SAMPLE"], cid)
        if key in emitted_complex:
            # second member of the pair will be appended when the first was emitted
            continue

        # Emit both rows of the complex at the first occurrence, preserving within-pair order
        pair = buckets[key]
        cat  = per_complex_category.get(key, "other")
        annotated = []
        for pr in pair:
            pr2 = dict(pr)
            pr2["CATEGORY"] = cat
            annotated.append(pr2)
        out_rows.extend(annotated)
        emitted_complex.add(key)

    # Write with counts header
    with args.output.open("w", newline="") as fh:
        fh.write("# category\tcount\n")
        for cat in ordered:
            fh.write(f"# {cat}\t{counts[cat]}\n")
        fh.write("#\n")
        wr = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        wr.writeheader()
        wr.writerows(out_rows)

    print(f"Wrote {len(out_rows)} data rows → {args.output}")
    print("Counts (INS+SNV complexes only):")
    for cat in ordered:
        print(f"  {cat:<34} {counts[cat]}")

if __name__ == "__main__":
    main()
