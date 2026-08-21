"""Behaviour when the truth itself is split into phase sets.

A ground truth assembled per haplotype block, rather than chromosome-wide, gives
no shared frame across its own block boundary: two variants in different truth
phase sets carry no relative-orientation claim, so no pair spanning that boundary
may be scored.  Every other fixture in the suite leaves the truth unphased into a
single block, which leaves these branches unexercised.
"""
from evaluate_phasing import (baseline_block_map, evaluate_genes, evaluate_new_junctions,
                              evaluate_phasing, evaluate_svs, get_overlapping_sites,
                              load_annotations, load_phasing)


def _snp(chrom, pos, gt, ps=1):
    return (chrom, pos, 'ACGT'[pos % 4], 'TGCA'[pos % 4], gt, ps)


def _ins(chrom, pos, gt, ps=1):
    return (chrom, pos, 'A', 'ATTT', gt, ps)


def _pair(write_vcf, target_records, truth_records):
    target = load_phasing(write_vcf('t.vcf', target_records))
    truth = load_phasing(write_vcf('g.vcf', truth_records))
    return get_overlapping_sites(target, truth), target, truth


def test_switches_are_not_counted_across_a_truth_block_boundary(write_vcf):
    """The intersection block splits when either file's phase set changes."""
    t = [_snp('chr1', p, '0|1', 1) for p in (10, 20, 30, 40)]
    g = [_snp('chr1', 10, '0|1', 100), _snp('chr1', 20, '0|1', 100),
         _snp('chr1', 30, '1|0', 200), _snp('chr1', 40, '1|0', 200)]
    pairs, target, truth = _pair(write_vcf, t, g)
    res = evaluate_phasing(pairs, target, truth)
    assert res.intersection_blocks == 2
    assert (res.assessed_pairs, res.switches) == (2, 0)   # the boundary pair is not scored


def test_sv_side_is_dropped_when_the_truth_splits_it(write_vcf):
    """The SV shares a target phase set with both anchors but not a truth one."""
    t = [_snp('chr1', 10, '0|1', 1), _ins('chr1', 20, '0|1', 1), _snp('chr1', 30, '0|1', 1)]
    g = [_snp('chr1', 10, '0|1', 100), _ins('chr1', 20, '0|1', 200),
         _snp('chr1', 30, '0|1', 200)]
    pairs, target, truth = _pair(write_vcf, t, g)
    res = evaluate_svs(pairs, target, truth)
    assert res.svs == 1
    assert res.one_sided == 1          # only the right anchor shares a truth block
    assert res.connections == 1
    assert res.correct + res.flipped + res.ambiguous == 0


def test_sv_with_both_sides_split_by_the_truth_has_no_anchor(write_vcf):
    t = [_snp('chr1', 10, '0|1', 1), _ins('chr1', 20, '0|1', 1), _snp('chr1', 30, '0|1', 1)]
    g = [_snp('chr1', 10, '0|1', 100), _ins('chr1', 20, '0|1', 200),
         _snp('chr1', 30, '0|1', 300)]
    pairs, target, truth = _pair(write_vcf, t, g)
    res = evaluate_svs(pairs, target, truth)
    assert (res.svs, res.no_anchor, res.connections) == (1, 1, 0)


def test_gene_pair_split_by_the_truth_is_unresolved(write_vcf, write_bed):
    t = [_snp('chr1', p, '0|1', 1) for p in (10, 20, 30)]
    g = [_snp('chr1', 10, '0|1', 100), _snp('chr1', 20, '0|1', 100),
         _snp('chr1', 30, '0|1', 200)]
    pairs, target, truth = _pair(write_vcf, t, g)
    exons = load_annotations(write_bed('e.bed', [('chr1', 0, 100, 'GENE1')]))
    res = evaluate_genes(pairs, target, truth, exons)
    assert (res.connections, res.unresolved, res.errors) == (1, 1, 0)
    assert res.genes_unresolved == 1
    assert res.connections + res.unresolved == res.sites - res.genes


def test_new_connection_split_by_the_truth_counts_as_no_truth_frame(write_vcf):
    """The join was made, but the truth cannot say whether it is right."""
    t = [_snp('chr1', p, '0|1', 1) for p in (10, 20, 30, 40)]
    b = [_snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0|1', 5),
         _snp('chr1', 30, '0|1', 6), _snp('chr1', 40, '0|1', 6)]
    g = [_snp('chr1', 10, '0|1', 100), _snp('chr1', 20, '0|1', 100),
         _snp('chr1', 30, '0|1', 200), _snp('chr1', 40, '0|1', 200)]
    target = load_phasing(write_vcf('t.vcf', t))
    baseline = load_phasing(write_vcf('b.vcf', b))
    truth = load_phasing(write_vcf('g.vcf', g))
    pairs = get_overlapping_sites(target, truth)
    res = evaluate_new_junctions(pairs, target, truth, baseline_block_map(target, baseline))

    assert res.structural_joins == 1
    assert res.junctions == 0            # nothing scorable
    assert res.no_truth_frame == 1       # and the reason is recorded
    shortfall = res.structural_joins - res.junctions - res.no_truth_frame
    assert shortfall == 0                # so it is not double counted as a missing variant


def test_chromosome_wide_also_merges_the_truth_blocks(write_vcf):
    """Documents a trap: `ignore_phase_blocks` concatenates blocks in BOTH files.

    The within-block statistic correctly declines to score the pair spanning the
    truth's own block boundary.  The chromosome-wide statistic merges the truth's
    blocks as well, which manufactures a frame that does not exist and reports the
    orientation difference across that boundary as a switch error.

    Harmless for a truth phased chromosome-wide -- the Q100 VCFs carry no PS at
    all, so every truth block is already one -- but it makes the chromosome-wide
    row unusable against a truth phased per assembly block.  Restricting the flag
    to the target would leave the numbers unchanged on a PS-free truth.
    """
    t = [_snp('chr1', p, '0|1', 1) for p in (10, 20, 30, 40)]
    g = [_snp('chr1', 10, '0|1', 100), _snp('chr1', 20, '0|1', 100),
         _snp('chr1', 30, '1|0', 200), _snp('chr1', 40, '1|0', 200)]
    pairs, target, truth = _pair(write_vcf, t, g)

    assert evaluate_phasing(pairs, target, truth).switches == 0
    assert evaluate_phasing(pairs, target, truth, ignore_phase_blocks=True).switches == 1

    # and with a truth that is genuinely one block, the two agree
    g_single = [_snp('chr1', 10, '0|1', 100), _snp('chr1', 20, '0|1', 100),
                _snp('chr1', 30, '1|0', 100), _snp('chr1', 40, '1|0', 100)]
    pairs2, target2, truth2 = _pair(write_vcf, t, g_single)
    assert evaluate_phasing(pairs2, target2, truth2).switches == 1
    assert evaluate_phasing(pairs2, target2, truth2, ignore_phase_blocks=True).switches == 1
