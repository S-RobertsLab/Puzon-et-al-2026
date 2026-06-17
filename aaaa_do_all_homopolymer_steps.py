#!/usr/bin/env python3
"""
run_homopolymer_pipeline.py
---------------------------
Pipeline driver that runs, in order:

  0) (built-in) complex-only prefilter       (keep rows with COMPLEX_ID != ".")
  1) filter_snv_ins_complex.py               (keep 2-row complexes with SNV + single-base INS)
  2) correct_snv_ins_complex.py              (right-justify INS before G-SNV and swap bases)
  3) flip_g_ins_snv.py                       (reverse-complement G-based INS+SNV to C-based on - strand)
  4) classify_homopolymers.py                (annotate categories and prepend counts)

Usage
-----
  python run_homopolymer_pipeline.py  input.complex
      [--outdir OUTDIR]
      [--keep-intermediate]
      [--scripts-dir DIR]
      [--prefix PREFIX]

Defaults
--------
- Outputs go next to input unless --outdir is specified.
- Intermediates are deleted unless --keep-intermediate is set.
- Scripts are searched in the same folder as this driver; override with --scripts-dir.

Final output:
  <outdir>/<prefix><stem>.snv_ins.fixed.flipped.classified.tsv
"""

from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run INS+SNV homopolymer pipeline (expects a .complex file).")
    p.add_argument("infile", type=Path, help="Input .complex TSV (with header)")
    p.add_argument("--outdir", type=Path, help="Output directory (default: input dir)")
    p.add_argument("--keep-intermediate", action="store_true",
                   help="Keep intermediate files (default: delete)")
    p.add_argument("--scripts-dir", type=Path,
                   help="Directory containing the component scripts "
                        "(default: same dir as this driver)")
    p.add_argument("--prefix", type=str, default="",
                   help="Optional prefix for output filenames")
    return p.parse_args()

def exe() -> str:
    """Return the current Python executable."""
    return sys.executable or "python3"

def ensure_exists(p: Path, what: str):
    if not p.exists():
        sys.exit(f"[FATAL] {what} not found: {p}")

def run(cmd: list[str], label: str):
    print(f"[RUN] {label}: {' '.join(map(str, cmd))}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"[FATAL] Step '{label}' failed with exit code {e.returncode}")

# ────────────────────────────────────────────────────────────────────────────
# Step 0: built-in complex-only prefilter (header-preserving, robust to column order)
# ────────────────────────────────────────────────────────────────────────────
def write_complex_only(infile: Path, outfile: Path):
    """
    Copy header + only rows where COMPLEX_ID != "."
    (Assumes tab-separated with a header line containing 'COMPLEX_ID'.)
    """
    with infile.open() as fi, outfile.open("w") as fo:
        header = fi.readline()
        if not header:
            sys.exit(f"[FATAL] Empty input: {infile}")
        fo.write(header)

        cols = header.rstrip("\n").split("\t")
        try:
            idx = cols.index("COMPLEX_ID")
        except ValueError:
            sys.exit("[FATAL] 'COMPLEX_ID' column not found in header.")

        for line in fi:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            # tolerate short/malformed lines
            if len(fields) <= idx:
                continue
            if fields[idx] != ".":
                fo.write(line)

def main():
    args = parse_args()

    # Resolve script locations
    driver_dir  = Path(__file__).parent.resolve()
    scripts_dir = (args.scripts_dir or driver_dir).resolve()

    # External scripts (use your filenames here)
    filter_script   = scripts_dir / 'filter_snv_ins_complex.py'
    correct_script  = scripts_dir / 'convert_miscalled_a3a_complex.py'
    flip_script     = scripts_dir / 'flip_fixed_SNV_INS_to_cbase.py'
    classify_script = scripts_dir / 'count_complex_snv_ins_homopolymer.py'

    for s, name in [(filter_script, "filter"),
                    (correct_script, "correct"),
                    (flip_script, "flip"),
                    (classify_script, "classify")]:
        ensure_exists(s, f"{name} script")

    # I/O paths
    infile = args.infile.resolve()
    ensure_exists(infile, "Input file (.complex)")
    outdir = (args.outdir or infile.parent).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    stem   = infile.stem  # if file is foo.complex, stem == 'foo'
    prefix = args.prefix

    # Intermediates (stable names)
    f0 = outdir / f"{prefix}{stem}.only.complex"                 # new prefilter output
    f1 = outdir / f"{prefix}{stem}.snv_ins.complex"
    f2 = outdir / f"{prefix}{stem}.snv_ins.fixed"
    f3 = outdir / f"{prefix}{stem}.snv_ins.fixed.flipped"
    f4 = outdir / f"{prefix}{stem}.snv_ins.fixed.flipped.classified.tsv"  # final

    # 0) Complex-only prefilter (header-preserving)
    print("[INFO] Step 0: complex-only prefilter")
    write_complex_only(infile, f0)

    # 1) Filter: keep 2-row complexes with SNV + single-base INS
    # run([exe(), str(filter_script), str(f0), "-o", str(f1)],
    #     "Filter (SNV + single-base INS)")

    run([exe(), str(filter_script), str(f0), "-o", str(f1),
     "--phasing", "assume-phased"],
    "Filter (SNV + single-base INS)")

    # 2) Correct: right-justify INS before G-SNV & swap bases when overlapping
    run([exe(), str(correct_script), str(f1), "-o", str(f2)],
        "Correct (justify & swap)")

    # 3) Flip: reverse-complement G-based complexes to C-based on minus strand
    run([exe(), str(flip_script), str(f2), "-o", str(f3)],
        "Flip to minus strand (C-based)")

    # 4) Classify homopolymers & prepend counts
    run([exe(), str(classify_script), str(f3), "-o", str(f4)],
        "Classify homopolymers")

    print("\n[OK] Pipeline complete.")
    print(f"     Final output: {f4}")

    # Cleanup
    if not args.keep_intermediate:
        for p in (f0, f1, f2, f3):
            try:
                p.unlink(missing_ok=True)
            except Exception:
                pass
    else:
        print("     (Intermediates kept)")

if __name__ == "__main__":
    main()
