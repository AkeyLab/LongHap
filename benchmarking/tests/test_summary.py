"""Formatting helpers and the summary row assembled by collect_summary."""
import argparse

import pytest

from evaluate_phasing import (collect_summary, evaluate_junctions, evaluate_phasing,
                              fraction, get_overlapping_sites, load_phasing, percent,
                              switch_columns, write_summary_tsv)


@pytest.mark.parametrize('num, den, expected', [
    (1, 4, '25.00%'),
    (0, 4, '0.00%'),
    (1, 3, '33.33%'),
    (1, 0, '--'),
    (0, 0, '--'),
])
def test_percent(num, den, expected):
    assert percent(num, den) == expected


@pytest.mark.parametrize('num, den, expected', [
    (1, 4, 0.25),
    (0, 4, 0.0),
    (1, 0, ''),
    (0, 0, ''),
])
def test_fraction(num, den, expected):
    assert fraction(num, den) == expected


def test_percent_and_fraction_agree():
    for num, den in ((3, 7), (218, 97159), (1, 2083)):
        assert percent(num, den) == f'{fraction(num, den) * 100:.2f}%'


def test_switch_columns_shape_is_the_same_with_and_without_a_result():
    class Fake:
        assessed_pairs = 10
        switches = 2
        switch_flips = (1, 1)
        intersection_blocks = 4
        blocks_no_switch = 3
        blocks_no_flip = 2
        blocks_clean = 1
        hamming = 3
        covered_variants = 11
        diff_genotypes = 0

    filled = switch_columns('within', Fake)
    empty = switch_columns('within', None)
    assert filled.keys() == empty.keys()
    assert all(v == '' for v in empty.values())
    assert filled['within_switch_rate'] == 0.2
    assert filled['within_switchflip_switches'] == 1
    assert filled['within_switchflip_rate'] == 0.2
    assert filled['within_blocks'] == 4
    assert filled['within_blocks_no_switch_rate'] == 0.75
    assert filled['within_blocks_clean_rate'] == 0.25


def _snp(chrom, pos, gt, ps=1):
    return (chrom, pos, 'ACGT'[pos % 4], 'TGCA'[pos % 4], gt, ps)


@pytest.fixture
def summary_inputs(write_vcf):
    target = load_phasing(write_vcf('t.vcf', [
        _snp('chr1', 10, '0|1'), _snp('chr1', 20, '0|1'), _snp('chr2', 30, '0|1')],
        contigs=('chr1', 'chr2')))
    truth = load_phasing(write_vcf('g.vcf', [
        _snp('chr1', 10, '0|1', None), _snp('chr1', 20, '1|0', None),
        _snp('chr2', 30, '0|1', None)], contigs=('chr1', 'chr2')))
    pairs = get_overlapping_sites(target, truth)
    args = argparse.Namespace(
        label='L', sample=None, vcf='t.vcf', gt_vcf='g.vcf', baseline_vcf=None,
        annotations=None, match='strict', only_snvs=False,
        keep_duplicate_positions=False, min_sv_length=0)
    return args, target, truth, pairs


def test_collect_summary_key_set_does_not_depend_on_optional_tests(summary_inputs):
    args, target, truth, pairs = summary_inputs
    within = evaluate_phasing(pairs, target, truth)
    chrwide = evaluate_phasing(pairs, target, truth, ignore_phase_blocks=True)
    junctions = evaluate_junctions(pairs, target, truth)

    row = collect_summary(args, target, truth, pairs, within, junctions, chrwide)
    assert row['label'] == 'L'
    assert row['chromosomes'] == 'chr1,chr2'
    assert row['within_switches'] == within.switches
    assert row['junction_assessed'] == junctions.junctions
    # optional sections still contribute their columns, empty
    for prefix in ('newconn_', 'sv_', 'gene_'):
        keys = [k for k in row if k.startswith(prefix)]
        assert keys and all(row[k] == '' for k in keys), prefix


def test_collect_summary_reports_every_chromosome_once(write_vcf):
    target = load_phasing(write_vcf('t.vcf', [
        _snp('chr1', 10, '0|1'), _snp('chr2', 20, '0|1'), _snp('chr1', 30, '0|1')],
        contigs=('chr1', 'chr2')))
    truth = load_phasing(write_vcf('g.vcf', [
        _snp('chr1', 10, '0|1', None)], contigs=('chr1', 'chr2')))
    pairs = get_overlapping_sites(target, truth)
    args = argparse.Namespace(
        label='', sample=None, vcf='t', gt_vcf='g', baseline_vcf=None, annotations=None,
        match='strict', only_snvs=False, keep_duplicate_positions=False, min_sv_length=0)
    within = evaluate_phasing(pairs, target, truth)
    row = collect_summary(args, target, truth, pairs, within,
                          evaluate_junctions(pairs, target, truth), within)
    assert row['chromosomes'] == 'chr1,chr2'      # in first-seen order, no repeats


def test_write_summary_tsv_round_trip(tmp_path):
    row = {'a': 1, 'b': '', 'c': 0.25}
    path = tmp_path / 'out.tsv'
    write_summary_tsv(str(path), row)
    header, values = (line.split('\t') for line in path.read_text().rstrip('\n').split('\n'))
    assert header == ['a', 'b', 'c']
    assert values == ['1', '', '0.25']


def test_write_summary_tsv_overwrites(tmp_path):
    path = tmp_path / 'out.tsv'
    write_summary_tsv(str(path), {'a': 1})
    write_summary_tsv(str(path), {'a': 2})
    assert path.read_text() == 'a\n2\n'
