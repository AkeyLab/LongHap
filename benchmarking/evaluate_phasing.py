#!/usr/bin/env python3
import argparse
import sys
from cyvcf2 import VCF
import numpy as np


def load_phasing(vcf_f):
    """
    Load phased heterozygous variants from a VCF file for a specific chromosome.
    """
    variant_type = []
    positions = []
    genotypes = []
    phased = []
    phase_block = []
    allele_a = []
    allele_b = []
    variant_calls = VCF(vcf_f)

    # only consider heterozygous variants
    for variant in variant_calls:
        if variant.gt_types == 1 and -1 not in variant.genotypes[0][:2]:
            pos = variant.POS
            ref = variant.REF.upper()
            alt = [a.upper() for a in variant.ALT]
            alleles = [ref] + alt
            gt = variant.genotypes[0][:2]
            al_a = alleles[gt[0]]
            al_b = alleles[gt[1]]
            is_phased = variant.genotypes[0][2]
            ps = variant.format('PS')[0][0] if 'PS' in variant.FORMAT else np.nan
            if len(alt) > 1:
                variant_type.append('Multi')

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
            allele_a.append(al_a)
            allele_b.append(al_b)

    variant_type = np.array(variant_type)
    positions = np.array(positions)
    genotypes = np.array(genotypes)
    phased = np.array(phased)
    phase_block = np.array(phase_block)
    allele_a = np.array(allele_a)
    allele_b = np.array(allele_b)
    return variant_type, positions, genotypes, phased, phase_block, allele_a, allele_b


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True, help='VCF to evaluate')
    parser.add_argument("--gt_vcf", required=True, help='Ground truth VCF')
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
     target_phased, target_phase_block, target_allele_a, target_allele_b) = load_phasing(args.vcf)
    (gt_variant_type, gt_positions, gt_genotypes,
     gt_phased, gt_phase_block, gt_allele_a, gt_allele_b) = load_phasing(args.gt_vcf)
    breakpoint()



if __name__ == '__main__':
    main(sys.argv[1:])

