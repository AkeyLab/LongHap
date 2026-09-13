"""The fixtures themselves.

If these fail, every other failure in the suite is suspect, so they come first.
"""
import pysam

from synthetic import ReadSpec, make_sequence, simulate_read, snv, Variant


def test_make_sequence_is_deterministic_and_low_homopolymer():
    a = make_sequence(500, seed=3)
    assert a == make_sequence(500, seed=3)
    assert make_sequence(500, seed=4) != a
    assert set(a) <= set('ACGT')
    for base in 'ACGT':
        assert base * 4 not in a


def test_simulate_read_snv():
    ref = 'A' * 20
    variants = [snv(11, 'A', 'G')]
    q0, c0 = simulate_read(ref, 0, 20, variants, hap=0)
    q1, c1 = simulate_read(ref, 0, 20, variants, hap=1)
    assert (q0, c0) == (ref, [(0, 20)])
    assert q1 == 'A' * 10 + 'G' + 'A' * 9
    assert c1 == [(0, 20)]


def test_simulate_read_insertion():
    ref = make_sequence(40, seed=1)
    var = Variant(11, ref[10], [ref[10] + 'TTT'])
    q, c = simulate_read(ref, 0, 40, [var], hap=1)
    assert c == [(0, 11), (1, 3), (0, 29)]
    assert q == ref[:11] + 'TTT' + ref[11:]
    assert len(q) == 43


def test_simulate_read_deletion():
    ref = make_sequence(40, seed=1)
    var = Variant(11, ref[10:14], [ref[10]])
    q, c = simulate_read(ref, 0, 40, [var], hap=1)
    assert c == [(0, 11), (2, 3), (0, 26)]
    assert q == ref[:11] + ref[14:]
    assert len(q) == 37


def test_simulate_read_skips_variants_outside_the_span():
    ref = make_sequence(100, seed=1)
    var = snv(51, ref[50], 'G' if ref[50] != 'G' else 'C')
    q, c = simulate_read(ref, 60, 90, [var], hap=1)
    assert (q, c) == (ref[60:90], [(0, 30)])


def test_write_reference_is_indexed(write_reference):
    seq = make_sequence(300, seed=2)
    path, written = write_reference(seq)
    assert written == seq
    with pysam.FastaFile(path) as fa:
        assert fa.fetch('chr1') == seq


def test_write_vcf_is_queryable(write_vcf):
    from cyvcf2 import VCF
    path = write_vcf([snv(100, 'A', 'G'), snv(200, 'C', 'T')], sequence_length=1000)
    got = [(v.POS, v.REF, v.ALT[0], v.gt_types[0]) for v in VCF(path)('chr1')]
    assert got == [(100, 'A', 'G', 1), (200, 'C', 'T', 1)]


def test_write_bam_round_trips(write_bam):
    ref = make_sequence(500, seed=5)
    var = snv(101, ref[100], 'G' if ref[100] != 'G' else 'C')
    path = write_bam(ref, [ReadSpec(0, 200, 0), ReadSpec(50, 250, 1)], [var])
    with pysam.AlignmentFile(path) as bam:
        reads = list(bam.fetch('chr1'))
    assert [r.reference_start for r in reads] == [0, 50]
    assert reads[0].query_sequence[100] == ref[100]
    assert reads[1].query_sequence[100 - 50] == var.alt[0]
    assert all(r.mapping_quality == 60 for r in reads)


def test_methylation_tags_round_trip(write_bam):
    ref = make_sequence(400, seed=6)
    c_positions = [i for i, b in enumerate(ref[:200]) if b == 'C'][:5]
    spec = ReadSpec(0, 200, 0, methylation={p: 0.9 for p in c_positions})
    path = write_bam(ref, [spec], [])
    with pysam.AlignmentFile(path) as bam:
        read = next(bam.fetch('chr1'))
    mods = read.modified_bases_forward
    assert ('C', 0, 'm') in mods
    called = dict(mods[('C', 0, 'm')])
    assert sorted(called) == c_positions
    assert all(p == 230 for p in called.values())


def test_make_locus_wires_everything(make_locus):
    ref_len = 2000
    variants = [snv(pos, 'A', 'G') for pos in (300, 500, 700)]
    locus = make_locus(variants, length=ref_len)
    with pysam.FastaFile(locus.ref_path) as fa:
        assert len(fa.fetch('chr1')) == ref_len
    # the VCF's REF must match the reference, so make_locus tests below build
    # their variants from the sequence; this one only checks plumbing exists.
    with pysam.AlignmentFile(locus.bam_path) as bam:
        assert bam.count('chr1') > 0
