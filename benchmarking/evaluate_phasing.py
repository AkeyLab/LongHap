#!/usr/bin/env python3
import argparse
import sys
from cyvcf2 import VCF
import numpy as np


def load_phasing(vcf_f, chrom):
    """
    Load phased heterozygous variants from a VCF file for a specific chromosome.
    """
    variant_type = []
    reference = []
    alternative = []
    positions = []
    genotypes = []
    phased = []
    phase_block = []

    variant_calls = VCF(vcf_f)

    # only consider heterozygous variants
    for variant in variant_calls(chrom):
        if variant.gt_types == 1 and -1 not in variant.genotypes[0][:2]:
            pos = variant.POS
            ref = variant.REF.upper()
            alt = [a.upper() for a in variant.ALT]
            gt = variant.genotypes[0][:2]
            is_phased = variant.genotypes[0][2]
            ps = variant.FORMAT['PS'] if 'PS' in variant.FORMAT else np.nan
            if len(alt) > 1:
                raise ValueError("Multi-allelic variants are not supported, "
                                 "please split the VCF into bi-allelic variants first, using bcftools norm -m -any")

            elif len(ref) == len(alt[0]) == 1:
                variant_type.append('SNP')
            elif len(ref) > len(alt[0]):
                variant_type.append('DEL')
            elif len(ref) < len(alt[0]):
                variant_type.append('INS')
            else:
                variant_type.append("Other")
            positions.append(pos)
            genotypes.append(gt)
            phased.append(is_phased)
            phase_block.append(ps)
            reference.append(ref)
            alternative.append(alt[0])

    variant_type = np.array(variant_type)
    positions = np.array(positions)
    genotypes = np.array(genotypes)
    phased = np.array(phased)
    phase_block = np.array(phase_block)
    reference = np.array(reference)
    alternative = np.array(alternative)
    return variant_type, positions, genotypes, phased, phase_block, reference, alternative


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, help='VCF to evaluate')
    parser.add_argument("--gt_vcf", required=True, help='Ground truth VCF')
    parser.add_argument("--chrom", required=True, help='Chromosome to evaluate')
    parser.add_argument('--evaluate_sv', action='store_true', default=False,
                        help='Evaluate structural variants instead of SNPs')
    parser.add_argument('--baseline_vcf', required=False,
                        help='Baseline VCF for evaluating new connection, for example, '
                             'derived from methylation information, in target VCF.')
    parser.add_argument('--annotations', required=False,
                        help='Bed file with genome annotations, e.g., CDS, '
                             'for which to evaluate the phasing of compound heterozygous variants')
    args = parser.parse_args()

    (target_variant_type, target_positions, target_genotypes,
     target_phased, target_phase_block, target_reference, target_alternative) = load_phasing(args.vcf, args.chrom)
    (gt_variant_type, gt_positions, gt_genotypes,
     gt_phased, gt_phase_block, gt_reference, gt_alternative) = load_phasing(args.gt_vcf, args.chrom)
    breakpoint()



if __name__ == '__main__':
    main(sys.argv[1:])

