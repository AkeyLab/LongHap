"""genes_by_variant and evaluate_genes: compound-heterozygosity scoring."""
import numpy as np
import pytest

from evaluate_phasing import (evaluate_genes, genes_by_variant, get_overlapping_sites,
                              load_annotations, load_phasing)


class FakeTarget:
    """Only what genes_by_variant reads: position and ref."""
    def __init__(self, spans):
        self.position = np.array([p for p, _ in spans])
        self.ref = np.array(['N' * n for _, n in spans], dtype=object)


# a SNP at 0-based 100 (span [100,101)) and a 5 bp deletion at 200 (span [200,205))
SPANS = [(100, 1), (200, 5)]


@pytest.mark.parametrize('label, exons, expected', [
    ('exon ends exactly at the SNP start', [(90, 100, 'A')], set()),
    ('exon starts exactly at the SNP', [(100, 105, 'B')], {'B'}),
    ('exon starts one base after the SNP', [(101, 110, 'C')], set()),
    ('exon ends one base after the SNP start', [(90, 101, 'D')], {'D'}),
    ('exon meets the end of the deletion span', [(204, 210, 'E')], {'E'}),
    ('exon just past the deletion span', [(205, 210, 'F')], set()),
    ('exon just before the deletion span', [(195, 200, 'G')], set()),
    ('exon spanning both variants', [(50, 300, 'H')], {'H'}),
])
def test_exon_overlap_is_half_open(label, exons, expected):
    got = set(genes_by_variant([0, 1], FakeTarget(SPANS), sorted(exons)))
    assert got == expected, label


def test_isoform_exons_record_a_variant_once():
    exons = sorted([(95, 105, 'H'), (98, 110, 'H'), (99, 102, 'H')])
    assert genes_by_variant([0, 1], FakeTarget(SPANS), exons)['H'] == [0]


def test_indel_reaching_into_an_exon_counts():
    """Overlap uses the full reference span, not just the start base."""
    assert genes_by_variant([1], FakeTarget(SPANS), [(203, 300, 'X')])['X'] == [1]


def _snp(chrom, pos, gt, ps=1):
    return (chrom, pos, 'ACGT'[pos % 4], 'TGCA'[pos % 4], gt, ps)


def _evaluate(write_vcf, write_bed, target_records, truth_records, bed_rows):
    target = load_phasing(write_vcf('t.vcf', target_records, contigs=('chr1', 'chr2')))
    truth = load_phasing(write_vcf('g.vcf', truth_records, contigs=('chr1', 'chr2')))
    exons = load_annotations(write_bed('e.bed', bed_rows))
    pairs = get_overlapping_sites(target, truth)
    return evaluate_genes(pairs, target, truth, exons)


def _check_identities(res):
    assert res.connections + res.unresolved == res.sites - res.genes
    assert res.genes_correct + res.genes_with_error + res.genes_unresolved == res.genes


def test_gene_with_all_sites_in_phase(write_vcf, write_bed):
    t = [_snp('chr1', p, '0|1') for p in (10, 20, 30)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30)]
    res = _evaluate(write_vcf, write_bed, t, g, [('chr1', 0, 100, 'GENE1')])
    assert (res.genes, res.sites, res.connections, res.errors) == (1, 3, 2, 0)
    assert res.genes_correct == 1
    _check_identities(res)


def test_gene_with_one_wrong_connection(write_vcf, write_bed):
    t = [_snp('chr1', 10, '0|1'), _snp('chr1', 20, '0|1'), _snp('chr1', 30, '1|0')]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30)]
    res = _evaluate(write_vcf, write_bed, t, g, [('chr1', 0, 100, 'GENE1')])
    assert (res.connections, res.errors) == (2, 1)
    assert (res.genes_correct, res.genes_with_error) == (0, 1)
    assert res.positions == [('chr1', 20, 30)]
    _check_identities(res)


def test_sites_split_across_phase_sets_are_unresolved(write_vcf, write_bed):
    t = [_snp('chr1', 10, '0|1', 1), _snp('chr1', 20, '0|1', 1), _snp('chr1', 30, '0|1', 2)]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30)]
    res = _evaluate(write_vcf, write_bed, t, g, [('chr1', 0, 100, 'GENE1')])
    assert (res.connections, res.unresolved, res.errors) == (1, 1, 0)
    assert res.genes_unresolved == 1
    _check_identities(res)


def test_counts_are_independent_of_phasing(write_vcf, write_bed):
    """genes/sites/single_site_genes must not move when a site is left unphased."""
    bed = [('chr1', 0, 100, 'GENE1')]
    g = [_snp('chr1', p, '0|1', None) for p in (10, 20, 30)]
    phased = [_snp('chr1', p, '0|1') for p in (10, 20, 30)]
    partly = [_snp('chr1', 10, '0|1'), _snp('chr1', 20, '0/1'), _snp('chr1', 30, '0|1')]

    a = _evaluate(write_vcf, write_bed, phased, g, bed)
    b = _evaluate(write_vcf, write_bed, partly, g, bed)
    assert (a.genes, a.sites, a.single_site_genes) == (b.genes, b.sites, b.single_site_genes)
    assert a.sites_scorable == 3 and b.sites_scorable == 2
    assert (a.connections, a.unresolved) == (2, 0)
    assert (b.connections, b.unresolved) == (0, 2)   # both pairs touch the unphased site
    _check_identities(b)


def test_single_site_gene_makes_no_connection(write_vcf, write_bed):
    t = [_snp('chr1', 10, '0|1'), _snp('chr1', 90, '0|1')]
    g = [_snp('chr1', 10, '0|1', None), _snp('chr1', 90, '0|1', None)]
    res = _evaluate(write_vcf, write_bed, t, g,
                    [('chr1', 0, 50, 'GENE1'), ('chr1', 50, 100, 'GENE2')])
    assert (res.genes, res.single_site_genes) == (0, 2)
    _check_identities(res)


def test_same_gene_name_on_two_chromosomes_is_kept_separate(write_vcf, write_bed):
    t = [_snp('chr1', 10, '0|1'), _snp('chr1', 20, '0|1'),
         _snp('chr2', 10, '0|1'), _snp('chr2', 20, '0|1')]
    g = [_snp('chr1', 10, '0|1', None), _snp('chr1', 20, '0|1', None),
         _snp('chr2', 10, '0|1', None), _snp('chr2', 20, '0|1', None)]
    res = _evaluate(write_vcf, write_bed, t, g,
                    [('chr1', 0, 100, 'PAR1'), ('chr2', 0, 100, 'PAR1')])
    assert res.genes == 2 and res.sites == 4 and res.connections == 2
    assert sorted(row[0] for row in res.per_gene) == ['chr1', 'chr2']
    _check_identities(res)


def test_annotation_chromosome_absent_from_the_vcf(write_vcf, write_bed):
    t = [_snp('chr1', 10, '0|1'), _snp('chr1', 20, '0|1')]
    g = [_snp('chr1', 10, '0|1', None), _snp('chr1', 20, '0|1', None)]
    res = _evaluate(write_vcf, write_bed, t, g,
                    [('chr1', 0, 100, 'GENE1'), ('chrX', 0, 100, 'GENEX')])
    assert res.genes == 1
    _check_identities(res)


def test_variant_outside_every_exon_is_ignored(write_vcf, write_bed):
    t = [_snp('chr1', 10, '0|1'), _snp('chr1', 20, '0|1'), _snp('chr1', 500, '1|0')]
    g = [_snp('chr1', 10, '0|1', None), _snp('chr1', 20, '0|1', None),
         _snp('chr1', 500, '0|1', None)]
    res = _evaluate(write_vcf, write_bed, t, g, [('chr1', 0, 100, 'GENE1')])
    assert (res.sites, res.errors) == (2, 0)     # the flipped site is intergenic
    _check_identities(res)


def test_no_exons_at_all(write_vcf, write_bed):
    t = [_snp('chr1', 10, '0|1'), _snp('chr1', 20, '0|1')]
    g = [_snp('chr1', 10, '0|1', None), _snp('chr1', 20, '0|1', None)]
    res = _evaluate(write_vcf, write_bed, t, g, [('chr2', 0, 100, 'GENE1')])
    assert (res.genes, res.sites, res.single_site_genes) == (0, 0, 0)
    _check_identities(res)
