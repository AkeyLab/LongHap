"""Switch statistics: the string helpers and evaluate_phasing."""
import pytest

from evaluate_phasing import (compute_switch_flips, evaluate_phasing, get_overlapping_sites,
                              hamming, load_phasing, switch_encoding)


def test_switch_encoding():
    assert switch_encoding('0001011') == '001110'
    assert switch_encoding('0') == ''
    assert switch_encoding('0000') == '000'


def test_hamming():
    assert hamming('ABCD', 'AXCY') == 2
    assert hamming('', '') == 0


@pytest.mark.parametrize('a, b, expected', [
    ('00011', '00100', (1, 0)),      # one switch
    ('00011', '00111', (0, 1)),      # a flip is two adjacent switches
    ('000', '001', (1, 0)),
    ('0101', '0101', (0, 0)),
    ('00000', '01000', (0, 1)),      # one variant out of place: two adjacent switches = a flip
    ('00000', '01110', (2, 0)),      # a run of three: two separate switch points, not a flip
])
def test_compute_switch_flips(a, b, expected):
    assert compute_switch_flips(a, b) == expected


def _evaluate(write_vcf, target_records, truth_records, **kwargs):
    target = load_phasing(write_vcf('t.vcf', target_records,
                                    contigs=('chr1', 'chr2')))
    truth = load_phasing(write_vcf('g.vcf', truth_records, contigs=('chr1', 'chr2')))
    pairs = get_overlapping_sites(target, truth)
    return evaluate_phasing(pairs, target, truth, **kwargs), target, truth


def _snp(chrom, pos, gt, ps=1):
    return (chrom, pos, 'ACGT'[pos % 4], 'TGCA'[pos % 4], gt, ps)


def test_no_switches_when_phasings_agree(write_vcf):
    t = [_snp('chr1', p, '0|1') for p in (10, 20, 30, 40)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    res, _, _ = _evaluate(write_vcf, t, g)
    assert (res.assessed_pairs, res.switches) == (3, 0)
    assert res.hamming == 0


def test_uniformly_inverted_phasing_is_not_a_switch(write_vcf):
    """Absolute orientation is arbitrary; only relative orientation matters."""
    t = [_snp('chr1', p, '1|0') for p in (10, 20, 30, 40)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    res, _, _ = _evaluate(write_vcf, t, g)
    assert (res.assessed_pairs, res.switches) == (3, 0)
    assert res.hamming == 0          # min(4 disagreeing, 0 agreeing)


def test_one_planted_switch(write_vcf):
    t = [_snp('chr1', 10, '0|1'), _snp('chr1', 20, '0|1'),
         _snp('chr1', 30, '1|0'), _snp('chr1', 40, '1|0')]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    res, _, _ = _evaluate(write_vcf, t, g)
    assert (res.assessed_pairs, res.switches) == (3, 1)
    assert res.switch_flips == (1, 0)
    assert res.hamming == 2


def test_single_flipped_variant_is_a_flip_not_two_switches(write_vcf):
    t = [_snp('chr1', 10, '0|1'), _snp('chr1', 20, '1|0'),
         _snp('chr1', 30, '0|1'), _snp('chr1', 40, '0|1')]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    res, _, _ = _evaluate(write_vcf, t, g)
    assert res.switches == 2                  # hamming of the switch encodings
    assert res.switch_flips == (0, 1)         # decomposed as one flip
    assert res.hamming == 1


def test_switches_are_not_counted_across_a_block_boundary(write_vcf):
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr1', 30, '1|0', 2), _snp('chr1', 40, '1|0', 2)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    res, _, _ = _evaluate(write_vcf, t, g)
    assert res.intersection_blocks == 2
    assert (res.assessed_pairs, res.switches) == (2, 0)   # 1 pair per block

    merged, _, _ = _evaluate(write_vcf, t, g, ignore_phase_blocks=True)
    assert (merged.assessed_pairs, merged.switches) == (3, 1)


def test_switches_are_not_counted_across_a_chromosome_boundary(write_vcf):
    """Even with one PS value shared, chr1 and chr2 are separate chains."""
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr2', 30, '1|0', 1), _snp('chr2', 40, '1|0', 1)]
    g = [_snp('chr1', 10, '0|1', None), _snp('chr1', 20, '0|1', None),
         _snp('chr2', 30, '0|1', None), _snp('chr2', 40, '0|1', None)]
    res, _, _ = _evaluate(write_vcf, t, g, ignore_phase_blocks=True)
    assert res.intersection_blocks == 2
    assert (res.assessed_pairs, res.switches) == (2, 0)


def test_singleton_blocks_are_skipped(write_vcf):
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 2)]
    g = [_snp('chr1', 10, '0|1', None), _snp('chr1', 20, '0|1', None)]
    res, _, _ = _evaluate(write_vcf, t, g)
    assert res.intersection_blocks == 0
    assert res.assessed_pairs == 0


def test_unphased_and_incomplete_variants_are_excluded(write_vcf):
    t = [_snp('chr1', 10, '0|1'), _snp('chr1', 20, '0/1'),
         _snp('chr1', 30, '.|1'), _snp('chr1', 40, '0|1')]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    res, _, _ = _evaluate(write_vcf, t, g)
    assert res.covered_variants == 2 and res.assessed_pairs == 1


def test_switch_positions_are_one_based_intervals(write_vcf):
    t = [_snp('chr1', 10, '0|1'), _snp('chr1', 20, '1|0')]
    g = [_snp('chr1', 10, '0|1', None), _snp('chr1', 20, '0|1', None)]
    res, _, _ = _evaluate(write_vcf, t, g)
    assert res.switch_positions == [('chr1', 10, 20)]


def test_agreement_is_anchored_on_haplotype_0(write_vcf):
    """Comparing the wrong haplotype must not pass unnoticed.

    For a site where the two files share an allele, comparing target haplotype 1
    against truth haplotype 0 yields the exact complement of the correct
    comparison, and a uniformly complemented agreement vector has the same switch
    count.  The two only diverge at a site whose genotype shares no allele with
    the truth's haplotype 0, so the chain below mixes both kinds.
    """
    t = [_snp('chr1', 10, '0|1'),
         ('chr1', 20, 'A', 'G,T', '0|1', 1),      # target carries A and G
         _snp('chr1', 30, '0|1')]
    g = [_snp('chr1', 10, '0|1', None),
         ('chr1', 20, 'A', 'G,T', '2|1', None),   # truth haplotype 0 carries T
         _snp('chr1', 30, '0|1', None)]
    res, _, _ = _evaluate(write_vcf, t, g)
    assert res.covered_variants == 3
    assert res.switches == 2                      # agree, disagree, agree
    assert res.diff_genotypes == 1                # the allele pairs differ
