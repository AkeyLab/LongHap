"""Methylation-assisted phasing.

The locus is built the way longhap's methylation step is meant to help.  Two
clusters of heterozygous variants sit far enough apart that no read covers a
variant from both, so read-backed phasing has to stop between them.  The reads
themselves do *overlap*, though, over a stretch carrying haplotype-specific
methylation: haplotype 0 reads are methylated there and haplotype 1 reads are
not.  That shared, differentially methylated stretch is the only thing that can
tie the right-hand reads to the left-hand haplotypes.

Note this is an overlap, not a span: a read reaching from one cluster's variants
into the other's would link them directly and the methylation step would have
nothing left to do.
"""
import numpy as np
import pytest

from synthetic import ReadSpec, alternating_snvs, cpg_positions

LEFT = [600, 700, 800]
RIGHT = [2400, 2500, 2600]
POSITIONS = LEFT + RIGHT
OVERLAP = (1400, 1600)     # where the two read sets meet


def build_reads(ref_seq, methylate=True):
    """Twelve reads over each cluster, meeting in the overlap window."""
    sites = cpg_positions(ref_seq, OVERLAP[0], OVERLAP[1], 60)
    specs = []
    for i in range(24):
        hap = i % 2
        start, end = (400, 1600) if i < 12 else (1400, 2800)
        meth = None
        if methylate:
            in_read = [p for p in sites if start <= p < end]
            prob = 0.95 if hap == 0 else 0.05
            meth = {p: prob for p in in_read}
        specs.append(ReadSpec(start, end, hap, name=f'r{i}', methylation=meth))
    return specs, sites


@pytest.fixture
def methylation_locus(make_locus):
    def _make(methylate=True):
        holder = {}

        def build(seq):
            holder['seq'] = seq
            return alternating_snvs(seq, POSITIONS, '010101')

        specs_holder = {}

        class _Specs(list):
            pass

        # make_locus needs the reference before reads can be rendered, so the
        # variant callback stashes the sequence and reads are built after.
        locus = make_locus(build, length=3200, read_specs=[])
        specs, sites = build_reads(holder['seq'], methylate)
        specs_holder['sites'] = sites
        return locus, specs, sites, holder['seq']
    return _make


def test_variant_data_alone_cannot_bridge_the_gap(make_locus, run_longhap, write_bam):
    """Baseline: with no methylation the two clusters stay separate blocks."""
    holder = {}

    def build(seq):
        holder['seq'] = seq
        return alternating_snvs(seq, POSITIONS, '010101')

    locus = make_locus(build, length=3200, read_specs=[])
    specs, _ = build_reads(holder['seq'], methylate=False)
    locus.bam_path = write_bam(holder['seq'], specs, locus.variants, name='m0.bam')

    lh = run_longhap(locus, write=False)
    assert len(lh.block_ends) >= 2, (
        f'expected two blocks without methylation, got {lh.block_ends}')


def test_methylation_step_joins_the_two_clusters(
        make_locus, run_longhap, write_bam, write_methylation_bed):
    holder = {}

    def build(seq):
        holder['seq'] = seq
        return alternating_snvs(seq, POSITIONS, '010101')

    locus = make_locus(build, length=3200, read_specs=[])
    specs, sites = build_reads(holder['seq'], methylate=True)
    locus.bam_path = write_bam(holder['seq'], specs, locus.variants, name='m1.bam')
    bed = write_methylation_bed([(p, 24, 50.0) for p in sites])

    lh = run_longhap(locus, write=False, methylation_calls_f=bed)
    assert len(lh.block_ends) == 1, (
        f'methylation should have joined the two clusters, got {lh.block_ends}')


def test_methylation_discriminates_the_two_haplotypes(
        make_locus, run_longhap, write_bam, write_methylation_bed, monkeypatch):
    """The likelihoods separate cleanly; only the threshold rejects them.

    This is the evidence behind ISSUE-15: the per-read scores for the correct
    and incorrect haplotype are far apart, yet both sit below log10(0.5), so the
    assignment rule discards a decision it had already made correctly.
    """
    from longhap import LongHap
    scores = []
    original = LongHap.calculate_probability_of_reads_belonging_to_haplotype_based_on_methylation

    def spy(*args, **kwargs):
        out = original(*args, **kwargs)
        scores.append(np.asarray(out).ravel())
        return out

    monkeypatch.setattr(
        LongHap,
        'calculate_probability_of_reads_belonging_to_haplotype_based_on_methylation',
        staticmethod(spy))

    holder = {}

    def build(seq):
        holder['seq'] = seq
        return alternating_snvs(seq, POSITIONS, '010101')

    locus = make_locus(build, length=3200, read_specs=[])
    specs, sites = build_reads(holder['seq'], methylate=True)
    locus.bam_path = write_bam(holder['seq'], specs, locus.variants, name='m1b.bam')
    bed = write_methylation_bed([(p, 24, 50.0) for p in sites])
    run_longhap(locus, write=False, methylation_calls_f=bed)

    assert len(scores) >= 2, 'the methylation step never scored any read'
    p_hap1, p_hap2 = scores[0], scores[1]
    margin = np.abs(p_hap1 - p_hap2)
    assert margin.min() > 10, f'haplotypes are not separated: margins {margin}'
    assert max(p_hap1.max(), p_hap2.max()) < np.log10(0.5), (
        'summed log-likelihoods sit below the single-site threshold the code '
        f'compares them to: best score {max(p_hap1.max(), p_hap2.max())}, '
        f'threshold {np.log10(0.5)}')


def test_methylation_leaves_variant_phasing_within_a_cluster_intact(
        make_locus, run_longhap, write_bam, write_methylation_bed):
    holder = {}

    def build(seq):
        holder['seq'] = seq
        return alternating_snvs(seq, POSITIONS, '010101')

    locus = make_locus(build, length=3200, read_specs=[])
    specs, sites = build_reads(holder['seq'], methylate=True)
    locus.bam_path = write_bam(holder['seq'], specs, locus.variants, name='m2.bam')
    bed = write_methylation_bed([(p, 24, 50.0) for p in sites])

    lh = run_longhap(locus, write=False, methylation_calls_f=bed)
    hap = lh.haplotypes[0]
    truth = [v.hap[0] for v in locus.variants]
    # within each cluster the read data alone settles the phase
    for cluster in (slice(0, 3), slice(3, 6)):
        got, want = hap[cluster], truth[cluster]
        flipped = [1 - w for w in want]
        assert list(got) == list(want) or list(got) == flipped, \
            f'{list(got)} matches neither {list(want)} nor {flipped}'


def test_differentially_methylated_sites_are_reported(
        make_locus, run_longhap, write_bam, write_methylation_bed, tmp_path):
    holder = {}

    def build(seq):
        holder['seq'] = seq
        return alternating_snvs(seq, POSITIONS, '010101')

    locus = make_locus(build, length=3200, read_specs=[])
    specs, sites = build_reads(holder['seq'], methylate=True)
    locus.bam_path = write_bam(holder['seq'], specs, locus.variants, name='m3.bam')
    bed = write_methylation_bed([(p, 24, 50.0) for p in sites])
    out = str(tmp_path / 'dms.tsv')

    lh = run_longhap(locus, methylation_calls_f=bed,
                     output_differentially_methylated_sites=out)
    header = open(out).readline().split()
    assert header[:3] == ['chrom', 'start', 'end'], header


def test_reads_without_mm_tags_are_skipped_by_the_methylation_step(
        make_locus, run_longhap, write_bam, write_methylation_bed):
    """Half the reads carry no MM tag; the run must still complete."""
    holder = {}

    def build(seq):
        holder['seq'] = seq
        return alternating_snvs(seq, POSITIONS, '010101')

    locus = make_locus(build, length=3200, read_specs=[])
    specs, sites = build_reads(holder['seq'], methylate=True)
    for spec in specs[::3]:
        spec.methylation = None
    locus.bam_path = write_bam(holder['seq'], specs, locus.variants, name='m4.bam')
    bed = write_methylation_bed([(p, 24, 50.0) for p in sites])

    lh = run_longhap(locus, write=False, methylation_calls_f=bed)
    assert lh.transition_matrix.shape[2] == len(POSITIONS) - 1


def test_methylation_with_no_informative_sites_is_a_no_op(
        make_locus, run_longhap, write_bam, write_methylation_bed):
    """A BED whose sites all fail the coverage/ratio filter must change nothing."""
    holder = {}

    def build(seq):
        holder['seq'] = seq
        return alternating_snvs(seq, POSITIONS, '010101')

    locus = make_locus(build, length=3200, read_specs=[])
    specs, sites = build_reads(holder['seq'], methylate=True)
    locus.bam_path = write_bam(holder['seq'], specs, locus.variants, name='m5.bam')

    baseline = run_longhap(locus, write=False).transition_matrix.copy()
    # coverage below 10 and ratio outside (20, 80): every site is filtered out
    bed = write_methylation_bed([(p, 4, 95.0) for p in sites])
    got = run_longhap(locus, write=False, methylation_calls_f=bed).transition_matrix
    assert np.allclose(got, baseline)
