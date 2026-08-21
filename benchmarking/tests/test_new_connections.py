"""baseline_block_map and evaluate_new_junctions: the joins a tool newly made."""
from evaluate_phasing import (baseline_block_map, evaluate_new_junctions,
                              get_overlapping_sites, load_phasing)


def _snp(chrom, pos, gt, ps=1):
    return (chrom, pos, 'ACGT'[pos % 4], 'TGCA'[pos % 4], gt, ps)


def _setup(write_vcf, target_records, baseline_records, truth_records):
    target = load_phasing(write_vcf('t.vcf', target_records, contigs=('chr1', 'chr2')))
    baseline = load_phasing(write_vcf('b.vcf', baseline_records, contigs=('chr1', 'chr2')))
    truth = load_phasing(write_vcf('g.vcf', truth_records, contigs=('chr1', 'chr2')))
    pairs = get_overlapping_sites(target, truth)
    return evaluate_new_junctions(pairs, target, truth,
                                  baseline_block_map(target, baseline))


def test_baseline_block_map_skips_unphased_baseline_variants(write_vcf):
    target = load_phasing(write_vcf('t.vcf', [
        _snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1)]))
    baseline = load_phasing(write_vcf('b.vcf', [
        _snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0/1', 5)]))
    assert baseline_block_map(target, baseline) == {0: ('chr1', 5)}


def test_one_new_join_correctly_oriented(write_vcf):
    """The target merges two baseline blocks and gets the orientation right."""
    t = [_snp('chr1', p, '0|1', 1) for p in (10, 20, 30, 40)]
    b = [_snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0|1', 5),
         _snp('chr1', 30, '0|1', 6), _snp('chr1', 40, '0|1', 6)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    res = _setup(write_vcf, t, b, g)
    assert (res.junctions, res.errors) == (1, 0)
    assert res.preexisting == 2                      # the two within-baseline-block pairs
    assert res.structural_joins == 1
    assert res.baseline_blocks - res.target_blocks == res.junctions


def test_one_new_join_wrongly_oriented(write_vcf):
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr1', 30, '1|0', 1), _snp('chr1', 40, '1|0', 1)]
    b = [_snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0|1', 5),
         _snp('chr1', 30, '1|0', 6), _snp('chr1', 40, '1|0', 6)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    res = _setup(write_vcf, t, b, g)
    assert (res.junctions, res.errors) == (1, 1)
    assert res.positions == [('chr1', 20, 30)]


def test_join_the_baseline_already_made_is_not_new(write_vcf):
    t = [_snp('chr1', p, '0|1', 1) for p in (10, 20, 30)]
    b = [_snp('chr1', p, '0|1', 5) for p in (10, 20, 30)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30)]
    res = _setup(write_vcf, t, b, g)
    assert (res.junctions, res.errors, res.preexisting) == (0, 0, 2)
    assert res.structural_joins == 0


def test_boundary_the_target_did_not_join_is_not_counted(write_vcf):
    """Different baseline blocks, but the target keeps them apart too."""
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr1', 30, '0|1', 2), _snp('chr1', 40, '0|1', 2)]
    b = [_snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0|1', 5),
         _snp('chr1', 30, '0|1', 6), _snp('chr1', 40, '0|1', 6)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    res = _setup(write_vcf, t, b, g)
    assert res.junctions == 0


def test_unphased_baseline_site_is_stepped_over(write_vcf):
    """The walk-outward rule: compare the nearest sites the baseline does phase.

    The middle site is unphased in the baseline, so the join between baseline
    blocks 5 and 6 is still scored, using the sites either side of it.
    """
    t = [_snp('chr1', p, '0|1', 1) for p in (10, 20, 30)]
    b = [_snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0/1', 5),
         _snp('chr1', 30, '0|1', 6)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30)]
    res = _setup(write_vcf, t, b, g)
    assert res.skipped_no_baseline_phase == 1
    assert (res.junctions, res.errors) == (1, 0)
    assert res.positions == []


def test_unphased_baseline_site_stepped_over_with_an_error(write_vcf):
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr1', 30, '1|0', 1)]
    b = [_snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0/1', 5),
         _snp('chr1', 30, '1|0', 6)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30)]
    res = _setup(write_vcf, t, b, g)
    assert (res.junctions, res.errors) == (1, 1)
    assert res.positions == [('chr1', 10, 30)]     # interval spans the skipped site


def test_baseline_identical_to_target_yields_no_new_connections(write_vcf):
    t = [_snp('chr1', p, '0|1', 1) for p in (10, 20, 30)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30)]
    res = _setup(write_vcf, t, t, g)
    assert (res.junctions, res.errors, res.structural_joins) == (0, 0, 0)
    assert res.baseline_blocks == res.target_blocks


def test_no_join_counted_across_a_chromosome_boundary(write_vcf):
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr2', 30, '0|1', 1), _snp('chr2', 40, '0|1', 1)]
    b = [_snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0|1', 5),
         _snp('chr2', 30, '0|1', 6), _snp('chr2', 40, '0|1', 6)]
    g = [_snp('chr1', 10, '0|1', None), _snp('chr1', 20, '0|1', None),
         _snp('chr2', 30, '0|1', None), _snp('chr2', 40, '0|1', None)]
    res = _setup(write_vcf, t, b, g)
    assert res.junctions == 0
    # the target block key carries the chromosome, so the two are distinct blocks
    assert res.target_blocks == 2 and res.baseline_blocks == 2


def test_join_with_no_truth_variant_is_reported_not_absorbed(write_vcf):
    """A baseline block the truth cannot see makes a join unscorable.

    Three baseline blocks are merged into one target block, so two joins were
    made, but the truth has no variant in the middle block.  The chain therefore
    jumps straight from block 5 to block 7 and only one join can be scored; the
    other must surface as a shortfall rather than vanish into the denominator.
    """
    t = [_snp('chr1', p, '0|1', 1) for p in (10, 20, 30, 40, 50, 60)]
    b = [_snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0|1', 5),
         _snp('chr1', 30, '0|1', 6), _snp('chr1', 40, '0|1', 6),
         _snp('chr1', 50, '0|1', 7), _snp('chr1', 60, '0|1', 7)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 50, 60)]   # middle block missing
    res = _setup(write_vcf, t, b, g)

    assert res.structural_joins == 2          # counted over everything the target phases
    assert res.junctions == 1                 # only one is scorable
    assert res.no_truth_frame == 0
    shortfall = res.structural_joins - res.junctions - res.no_truth_frame
    assert shortfall == 1
    assert res.baseline_blocks - res.target_blocks == res.junctions


# --------------------------------------------------------------------------
# call sets that strand a variant at the old phase set
#
# methphaser relabels and flips a block it merges but leaves the block's first
# variant behind, still carrying the old PS and the un-flipped genotype.  The two
# variants either side of the join are then not adjacent on the chromosome, and
# the stranded one is a phase set of its own.
# --------------------------------------------------------------------------

def _merge_with_orphan(flip):
    """Baseline blocks 5 and 6 merged into target block 1, stranding site 40.

    The truth has block 6 inverted relative to block 5, so ``flip=True`` is the
    correct merge.  Either way site 40 keeps PS 6 and the genotype the baseline
    gave it, which is what the scan has to step over.
    """
    moved = '1|0' if flip else '0|1'
    target = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
              _snp('chr1', 30, '0|1', 1),
              _snp('chr1', 40, '0|1', 6),              # stranded at the old PS
              _snp('chr1', 50, moved, 1), _snp('chr1', 60, moved, 1)]
    baseline = [_snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0|1', 5),
                _snp('chr1', 30, '0|1', 5),
                _snp('chr1', 40, '0|1', 6), _snp('chr1', 50, '0|1', 6),
                _snp('chr1', 60, '0|1', 6)]
    truth = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30)] + \
            [_snp('chr1', p, '1|0', None) for p in (40, 50, 60)]
    return target, baseline, truth


def test_join_is_scored_although_a_variant_was_left_at_the_old_phase_set(write_vcf):
    """The join the target really made, across the stranded variant."""
    res = _setup(write_vcf, *_merge_with_orphan(flip=True))
    assert (res.junctions, res.errors) == (1, 0)
    assert res.structural_joins == 1


def test_a_wrong_join_across_a_stranded_variant_is_still_an_error(write_vcf):
    """And the verdict comes from the relabelled sites, not the stranded one.

    The stranded site 40 carries the un-flipped genotype, so a scan that scored
    against it would judge the merge on the baseline's arbitrary cross-block
    orientation.  The reported interval names the two sites that carry the
    decision, which is how the two cases are told apart.
    """
    res = _setup(write_vcf, *_merge_with_orphan(flip=False))
    assert (res.junctions, res.errors) == (1, 1)
    assert res.positions == [('chr1', 30, 50)]


def test_stranding_a_variant_changes_the_verdict_for_nobody(write_vcf):
    """The same merge with the PS written out properly is scored the same way.

    Only the reported interval moves, because it now stops at site 40 instead of
    spanning it -- the same widening the walk already does for a site the baseline
    leaves unphased.
    """
    stranded = _setup(write_vcf, *_merge_with_orphan(flip=False))
    target, baseline, truth = _merge_with_orphan(flip=False)
    target[3] = _snp('chr1', 40, '0|1', 1)             # site 40 relabelled as well
    intact = _setup(write_vcf, target, baseline, truth)
    assert (intact.junctions, intact.errors) == (stranded.junctions, stranded.errors)
    assert intact.structural_joins == stranded.structural_joins
    assert intact.positions == [('chr1', 30, 40)]
    assert stranded.positions == [('chr1', 30, 50)]
    assert intact.singleton_target_blocks == 0
    assert intact.baseline_blocks - intact.target_blocks == intact.junctions


def test_every_structural_join_is_scored_when_the_truth_can_frame_it(write_vcf):
    """The property the chromosome-wide walk violated: nothing silently dropped."""
    for flip in (True, False):
        res = _setup(write_vcf, *_merge_with_orphan(flip))
        assert res.junctions == res.structural_joins
        assert res.no_truth_frame == 0


def test_three_baseline_blocks_merged_with_two_stranded_variants(write_vcf):
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr1', 30, '0|1', 6),                    # stranded
         _snp('chr1', 40, '0|1', 1),
         _snp('chr1', 50, '0|1', 7),                    # stranded
         _snp('chr1', 60, '0|1', 1)]
    b = [_snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0|1', 5),
         _snp('chr1', 30, '0|1', 6), _snp('chr1', 40, '0|1', 6),
         _snp('chr1', 50, '0|1', 7), _snp('chr1', 60, '0|1', 7)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40, 50, 60)]
    res = _setup(write_vcf, t, b, g)
    assert (res.junctions, res.errors) == (2, 0)
    assert res.structural_joins == 2
    assert res.singleton_target_blocks == 2


def test_singleton_target_blocks_account_for_the_block_count_shortfall(write_vcf):
    """Why baseline_blocks - target_blocks stops equalling the joins made.

    Each stranded variant leaves a target block behind that cancels the one the
    merge removed, so the block-count difference under-reports the joins by
    exactly the number of blocks stranded that way.
    """
    res = _setup(write_vcf, *_merge_with_orphan(flip=True))
    assert res.singleton_target_blocks == 1
    assert res.baseline_blocks - res.target_blocks == 0        # not 1
    assert res.baseline_blocks - res.target_blocks + res.singleton_target_blocks \
        == res.junctions


def test_a_stranded_variant_the_baseline_leaves_unphased_creates_no_junction(write_vcf):
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1),
         _snp('chr1', 30, '0|1', 6),                    # stranded
         _snp('chr1', 40, '0|1', 1), _snp('chr1', 50, '0|1', 1)]
    b = [_snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0|1', 5),
         _snp('chr1', 30, '0/1', None),                 # and unphased in the baseline
         _snp('chr1', 40, '0|1', 6), _snp('chr1', 50, '0|1', 6)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40, 50)]
    res = _setup(write_vcf, t, b, g)
    assert res.skipped_no_baseline_phase == 1
    assert (res.junctions, res.errors) == (1, 0)
    assert res.positions == []
    assert res.singleton_target_blocks == 0       # it never entered the scan


# --------------------------------------------------------------------------
# the interleaving guard
# --------------------------------------------------------------------------

def test_interleaved_baseline_blocks_are_flagged(write_vcf):
    """Two baseline blocks alternating along one target block.

    Three baseline-block changes where merging two blocks can only be one join,
    so the scan would over-count.  Not seen in this data; the counter is here so
    that it cannot start happening unnoticed.
    """
    t = [_snp('chr1', p, '0|1', 1) for p in (10, 20, 30, 40)]
    b = [_snp('chr1', 10, '0|1', 5), _snp('chr1', 20, '0|1', 6),
         _snp('chr1', 30, '0|1', 5), _snp('chr1', 40, '0|1', 6)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30, 40)]
    res = _setup(write_vcf, t, b, g)
    assert res.interleaved_blocks == 1
    assert res.junctions == 3 > res.structural_joins == 1


def test_ordinary_data_is_not_flagged_as_interleaved(write_vcf):
    for records in (_merge_with_orphan(True), _merge_with_orphan(False)):
        assert _setup(write_vcf, *records).interleaved_blocks == 0
