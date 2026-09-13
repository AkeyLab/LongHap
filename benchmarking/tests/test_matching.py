"""get_overlapping_sites: how variants are paired across the two files."""
from evaluate_phasing import get_overlapping_sites, load_phasing


def test_strict_requires_ref_and_alt_to_agree(write_vcf):
    target = load_phasing(write_vcf('t.vcf', [
        ('chr1', 10, 'A', 'G', '0|1', 1),
        ('chr1', 20, 'C', 'T', '0|1', 1),
    ]))
    truth = load_phasing(write_vcf('g.vcf', [
        ('chr1', 10, 'A', 'G', '0|1', None),
        ('chr1', 20, 'C', 'A', '0|1', None),      # same position, different ALT
    ]))
    pairs = get_overlapping_sites(target, truth)
    assert pairs.tolist() == [[0, 0]]


def test_alleles_mode_pairs_an_insertion_with_a_deletion(write_vcf):
    """The chr13 case: C->CA in one file against CA->C in the other.

    The unordered allele pair is identical, so `alleles` matches them, but REF
    differs so they are not the same variant and `strict` rejects them.
    """
    target = load_phasing(write_vcf('t.vcf', [('chr1', 10, 'CA', 'C', '0|1', 1)]))
    truth = load_phasing(write_vcf('g.vcf', [('chr1', 10, 'C', 'CA', '1|0', None)]))

    assert get_overlapping_sites(target, truth, match='strict').tolist() == []
    assert get_overlapping_sites(target, truth, match='alleles').tolist() == [[0, 0]]


def test_alleles_mode_still_requires_the_same_allele_pair(write_vcf):
    target = load_phasing(write_vcf('t.vcf', [('chr1', 10, 'A', 'G', '0|1', 1)]))
    truth = load_phasing(write_vcf('g.vcf', [('chr1', 10, 'A', 'T', '0|1', None)]))
    assert get_overlapping_sites(target, truth, match='alleles').tolist() == []


def test_alleles_mode_skips_incomplete_genotypes(write_vcf):
    target = load_phasing(write_vcf('t.vcf', [('chr1', 10, 'A', 'G', '0|1', 1)]))
    truth = load_phasing(write_vcf('g.vcf', [('chr1', 10, 'A', 'G', '.|1', None)]))
    assert get_overlapping_sites(target, truth, match='alleles').tolist() == []
    # strict keys on the record, not the genotype, so it still pairs them
    assert get_overlapping_sites(target, truth, match='strict').tolist() == [[0, 0]]


def test_match_is_per_chromosome(write_vcf):
    target = load_phasing(write_vcf('t.vcf', [
        ('chr1', 10, 'A', 'G', '0|1', 1),
        ('chr2', 10, 'A', 'G', '0|1', 1),
    ], contigs=('chr1', 'chr2')))
    truth = load_phasing(write_vcf('g.vcf', [
        ('chr2', 10, 'A', 'G', '0|1', None),
    ], contigs=('chr1', 'chr2')))
    assert get_overlapping_sites(target, truth).tolist() == [[1, 0]]


def test_no_common_variants(write_vcf):
    target = load_phasing(write_vcf('t.vcf', [('chr1', 10, 'A', 'G', '0|1', 1)]))
    truth = load_phasing(write_vcf('g.vcf', [('chr1', 99, 'A', 'G', '0|1', None)]))
    pairs = get_overlapping_sites(target, truth)
    assert pairs.shape == (0, 2)      # empty but still 2-D, so indexing works
