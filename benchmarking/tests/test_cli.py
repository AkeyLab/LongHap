"""End-to-end runs of the script, and the summary TSV contract."""
import re

import pytest


def _snp(chrom, pos, gt, ps=1):
    return (chrom, pos, 'ACGT'[pos % 4], 'TGCA'[pos % 4], gt, ps)


def _ins(chrom, pos, gt, ps=1):
    return (chrom, pos, 'A', 'ATTT', gt, ps)


@pytest.fixture
def dataset(write_vcf, write_bed):
    """A small two-chromosome dataset exercising every metric."""
    target = write_vcf('t.vcf', [
        _snp('chr1', 10, '0|1', 1), _ins('chr1', 20, '0|1', 1),
        _snp('chr1', 30, '0|1', 1), _snp('chr1', 40, '1|0', 1),
        _snp('chr1', 50, '0|1', 2), _snp('chr1', 60, '0|1', 2),
        _snp('chr2', 10, '0|1', 3), _snp('chr2', 20, '0|1', 3),
    ], contigs=('chr1', 'chr2'))
    baseline = write_vcf('b.vcf', [
        _snp('chr1', 10, '0|1', 7), _ins('chr1', 20, '0|1', 7),
        _snp('chr1', 30, '0|1', 8), _snp('chr1', 40, '1|0', 8),
        _snp('chr1', 50, '0|1', 9), _snp('chr1', 60, '0|1', 9),
        _snp('chr2', 10, '0|1', 10), _snp('chr2', 20, '0|1', 10),
    ], contigs=('chr1', 'chr2'))
    truth = write_vcf('g.vcf', [
        _snp('chr1', 10, '0|1', None), _ins('chr1', 20, '0|1', None),
        _snp('chr1', 30, '0|1', None), _snp('chr1', 40, '0|1', None),
        _snp('chr1', 50, '0|1', None), _snp('chr1', 60, '0|1', None),
        _snp('chr2', 10, '0|1', None), _snp('chr2', 20, '0|1', None),
    ], contigs=('chr1', 'chr2'))
    exons = write_bed('e.bed', [('chr1', 0, 45, 'GENE1'), ('chr2', 0, 100, 'GENE2')])
    return target, baseline, truth, exons


def test_minimal_run(run_cli, dataset):
    target, _, truth, _ = dataset
    out = run_cli('--vcf', target, '--gt_vcf', truth).stdout
    assert 'WITHIN PHASE BLOCKS' in out
    assert 'CHROMOSOME-WIDE' in out
    assert 'SV PHASING' not in out            # optional sections stay off
    assert 'NEW CONNECTIONS' not in out
    assert 'COMPOUND HETEROZYGOSITY' not in out


def test_every_output_file_is_written(run_cli, dataset, tmp_path):
    target, baseline, truth, exons = dataset
    paths = {name: tmp_path / name for name in (
        'summary.tsv', 'genes.tsv', 'switch.bed', 'junction.bed', 'new.bed', 'sv.bed')}
    run_cli('--vcf', target, '--gt_vcf', truth, '--baseline_vcf', baseline,
            '--annotations', exons, '--evaluate_sv',
            '--summary_tsv', paths['summary.tsv'], '--gene_tsv', paths['genes.tsv'],
            '--switch_error_bed', paths['switch.bed'],
            '--junction_error_bed', paths['junction.bed'],
            '--new_connection_bed', paths['new.bed'], '--sv_error_bed', paths['sv.bed'])
    for name, path in paths.items():
        assert path.exists(), name
    assert paths['genes.tsv'].read_text().startswith('chrom\tgene\t')


def test_summary_tsv_is_two_lines_with_matching_widths(run_cli, dataset, tmp_path):
    target, baseline, truth, exons = dataset
    out = tmp_path / 'summary.tsv'
    run_cli('--vcf', target, '--gt_vcf', truth, '--baseline_vcf', baseline,
            '--annotations', exons, '--evaluate_sv', '--label', 'run1', '--summary_tsv', out)
    lines = out.read_text().rstrip('\n').split('\n')
    assert len(lines) == 2
    header, values = (line.split('\t') for line in lines)
    assert len(header) == len(values)
    row = dict(zip(header, values))
    assert row['label'] == 'run1'
    assert row['chromosomes'] == 'chr1,chr2'


def test_summary_header_is_stable_without_optional_tests(run_cli, dataset, tmp_path):
    target, baseline, truth, exons = dataset
    full, minimal = tmp_path / 'full.tsv', tmp_path / 'min.tsv'
    run_cli('--vcf', target, '--gt_vcf', truth, '--baseline_vcf', baseline,
            '--annotations', exons, '--evaluate_sv', '--summary_tsv', full)
    run_cli('--vcf', target, '--gt_vcf', truth, '--summary_tsv', minimal)

    full_header = full.read_text().split('\n')[0]
    assert full_header == minimal.read_text().split('\n')[0]

    row = dict(zip(full_header.split('\t'), minimal.read_text().split('\n')[1].split('\t')))
    assert row['sv_total'] == '' and row['gene_genes'] == '' and row['newconn_assessed'] == ''
    assert row['within_assessed_pairs'] != ''


def test_summary_tsv_overwrites_rather_than_appends(run_cli, dataset, tmp_path):
    target, _, truth, _ = dataset
    out = tmp_path / 'summary.tsv'
    run_cli('--vcf', target, '--gt_vcf', truth, '--summary_tsv', out)
    run_cli('--vcf', target, '--gt_vcf', truth, '--summary_tsv', out)
    assert len(out.read_text().rstrip('\n').split('\n')) == 2


def test_rates_are_fractions_or_empty(run_cli, dataset, tmp_path):
    target, baseline, truth, exons = dataset
    out = tmp_path / 'summary.tsv'
    run_cli('--vcf', target, '--gt_vcf', truth, '--baseline_vcf', baseline,
            '--annotations', exons, '--evaluate_sv', '--summary_tsv', out)
    header, values = (line.split('\t') for line in out.read_text().rstrip('\n').split('\n'))
    rates = {k: v for k, v in zip(header, values) if k.endswith('_rate')}
    assert rates
    for name, value in rates.items():
        if value == '':
            continue
        assert '%' not in value, name
        assert 0.0 <= float(value) <= 1.0, name


def test_every_reported_number_reaches_the_summary_tsv(run_cli, dataset, tmp_path):
    """The completeness contract: nothing printed may be missing from the TSV."""
    target, baseline, truth, exons = dataset
    out = tmp_path / 'summary.tsv'
    stdout = run_cli('--vcf', target, '--gt_vcf', truth, '--baseline_vcf', baseline,
                     '--annotations', exons, '--evaluate_sv', '--summary_tsv', out).stdout
    header, values = (line.split('\t') for line in out.read_text().rstrip('\n').split('\n'))
    row = dict(zip(header, values))

    printed = set()
    for line in stdout.splitlines():
        _, _, rhs = line.partition(':')
        for token in rhs.split():
            if re.fullmatch(r'-?\d+', token):
                printed.add(int(token))
            elif re.fullmatch(r'\d+/\d+', token):
                printed.update(int(x) for x in token.split('/'))
    in_tsv = {int(v) for v in row.values() if re.fullmatch(r'-?\d+', str(v))}
    assert printed <= in_tsv, f'missing from the summary TSV: {sorted(printed - in_tsv)}'

    # and every printed percentage must correspond to a rate column
    for line in stdout.splitlines():
        label, _, rhs = line.partition(':')
        for token in rhs.split():
            if re.fullmatch(r'\d+\.\d+%', token):
                pct = float(token[:-1])
                assert any(abs(float(v) * 100 - pct) < 0.005
                           for k, v in row.items() if k.endswith('_rate') and v != ''), label


def test_only_snvs_combines_with_the_other_tests(run_cli, dataset, tmp_path):
    target, baseline, truth, exons = dataset
    out = tmp_path / 'summary.tsv'
    stdout = run_cli('--vcf', target, '--gt_vcf', truth, '--baseline_vcf', baseline,
                     '--annotations', exons, '--evaluate_sv', '--only_snvs',
                     '--summary_tsv', out).stdout
    header, values = (line.split('\t') for line in out.read_text().rstrip('\n').split('\n'))
    row = dict(zip(header, values))
    assert row['only_snvs'] == '1'
    assert row['sv_total'] == '0'        # the indel was filtered out before typing
    assert 'SV PHASING' in stdout


def test_match_alleles_mode_runs(run_cli, dataset, tmp_path):
    target, _, truth, _ = dataset
    out = tmp_path / 'summary.tsv'
    run_cli('--vcf', target, '--gt_vcf', truth, '--match', 'alleles', '--summary_tsv', out)
    header, values = (line.split('\t') for line in out.read_text().rstrip('\n').split('\n'))
    assert dict(zip(header, values))['match'] == 'alleles'


def test_empty_input_does_not_crash(run_cli, write_vcf, tmp_path):
    empty = write_vcf('empty.vcf', [])
    other = write_vcf('other.vcf', [_snp('chr1', 10, '0|1')])
    out = tmp_path / 'summary.tsv'
    stdout = run_cli('--vcf', empty, '--gt_vcf', other, '--summary_tsv', out).stdout
    assert 'VARIANT COUNTS' in stdout
    assert len(out.read_text().rstrip('\n').split('\n')) == 2


def test_no_shared_chromosome_does_not_crash(run_cli, write_vcf, tmp_path):
    a = write_vcf('a.vcf', [_snp('chr1', 10, '0|1')], contigs=('chr1', 'chr2'))
    b = write_vcf('b.vcf', [_snp('chr2', 10, '0|1', None)], contigs=('chr1', 'chr2'))
    out = tmp_path / 'summary.tsv'
    run_cli('--vcf', a, '--gt_vcf', b, '--summary_tsv', out)
    header, values = (line.split('\t') for line in out.read_text().rstrip('\n').split('\n'))
    assert dict(zip(header, values))['common_het'] == '0'
