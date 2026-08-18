#!/usr/bin/env python
"""
Fix duplicate/contradictory records in methphaser (v0.0.3) output VCFs.

methphaser's post-processing step (get_altered_vcf in meth_phaser_post_processing)
mixes 1-based-inclusive whatshap block coordinates into pysam's 0-based
VariantFile.fetch() calls. The resulting off-by-one causes a block's first
variant to sometimes be written twice: once (correctly) under the merged
block's phase set, and once (incorrectly) as an orphaned singleton whose PS
tag equals its own position, carrying unflipped/unmerged phase.

Fix: for any (CHROM, POS, REF, ALT) that appears more than once, keep only
the record whose PS tag belongs to the largest phase set in the file (the
real, established merge decision) and drop the rest. Records with a
different REF/ALT at the same POS (legitimate overlapping variants, e.g.
an SNV and a nearby SV call) are left untouched.

Note: methphaser also silently drops variants at chromosome/block-list
boundaries (its unphased-gap fill list is built from the snp_phased_block_1
column shifted by one row, so it can cover neither the region before the
first phase block, nor the final inter-block gap, nor the region after the
last phase block). This script does NOT attempt to recover those records.

Usage:
    python dedup_methphaser_vcf.py <input.vcf.gz> <output.vcf.gz>

Requires: cyvcf2, bcftools (for indexing the output).
"""

import sys
import argparse
import subprocess
from collections import defaultdict

from cyvcf2 import VCF, Writer


def dedup_methphaser_vcf(in_path, out_path):
    """Drop spurious duplicate records, keeping the one from the larger PS group."""
    # pass 1: count how many variants share each (CHROM, PS) group
    vcf = VCF(in_path)
    ps_counts = defaultdict(int)
    n_total = 0
    for rec in vcf:
        n_total += 1
        ps = rec.format("PS")
        if ps is not None:
            ps_counts[(rec.CHROM, int(ps[0][0]))] += 1
    vcf.close()

    def group_size(rec):
        ps = rec.format("PS")
        if ps is None:
            return 0
        return ps_counts.get((rec.CHROM, int(ps[0][0])), 0)

    # pass 2: for each run of records sharing (CHROM, POS, REF, ALT),
    # keep only the one belonging to the largest PS group
    vcf = VCF(in_path)
    writer = Writer(out_path, vcf)

    n_dropped = 0

    def flush(buffer):
        nonlocal n_dropped
        if not buffer:
            return
        best = buffer[0] if len(buffer) == 1 else max(buffer, key=group_size)
        writer.write_record(best)
        n_dropped += len(buffer) - 1

    buffer, last_key = [], None
    for rec in vcf:
        key = (rec.CHROM, rec.POS, rec.REF, tuple(rec.ALT))
        if key == last_key:
            buffer.append(rec)
        else:
            flush(buffer)
            buffer, last_key = [rec], key
    flush(buffer)

    writer.close()
    vcf.close()

    return n_total, n_dropped


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main():
    parser = argparse.ArgumentParser(
        description="Drop duplicate/contradictory records from a methphaser VCF."
    )
    parser.add_argument("input_vcf", help="input VCF (bgzipped, .vcf.gz)")
    parser.add_argument("output_vcf", help="output VCF (bgzipped, .vcf.gz)")
    args = parser.parse_args()

    n_total, n_dropped = dedup_methphaser_vcf(args.input_vcf, args.output_vcf)
    run(["bcftools", "index", "-t", args.output_vcf])

    print(f"total input records: {n_total}", file=sys.stderr)
    print(f"duplicate records dropped: {n_dropped}", file=sys.stderr)
    print(f"output indexed and ready: {args.output_vcf}", file=sys.stderr)


if __name__ == "__main__":
    main()
