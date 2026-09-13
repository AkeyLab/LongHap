"""The command line entry point and the intermediate-file cache."""
import os
import subprocess
import sys
from pathlib import Path

import pytest

from synthetic import REPO_ROOT, alternating_snvs

SCRIPT = REPO_ROOT / 'longhap.py'
POSITIONS = list(range(400, 3600, 200))


@pytest.fixture
def cli_locus(make_locus):
    return make_locus(lambda seq: alternating_snvs(seq, POSITIONS))


def run_cli(locus, tmp_path, *extra, expect_ok=True):
    out = tmp_path / 'cli.vcf.gz'
    cmd = [sys.executable, str(SCRIPT),
           '--vcf', locus.vcf_path, '-b', locus.bam_path,
           '-r', locus.ref_path, '-c', locus.chrom,
           '-o', str(out), '--pacbio',
           '--log', str(tmp_path / 'longhap.log'), *map(str, extra)]
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(tmp_path))
    if expect_ok:
        assert proc.returncode == 0, f'exit {proc.returncode}\n{proc.stdout}\n{proc.stderr}'
    return proc, out


def test_cli_phases_and_writes_a_vcf(cli_locus, tmp_path, read_phased_vcf):
    _, out = run_cli(cli_locus, tmp_path)
    records = read_phased_vcf(str(out))
    assert [pos for pos, _, _ in records] == POSITIONS
    assert all('|' in gt for _, gt, _ in records)


def test_cli_requires_a_sequencing_technology(cli_locus, tmp_path):
    out = tmp_path / 'cli.vcf.gz'
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), '--vcf', cli_locus.vcf_path,
         '-b', cli_locus.bam_path, '-r', cli_locus.ref_path, '-c', cli_locus.chrom,
         '-o', str(out), '--log', str(tmp_path / 'l.log')],
        capture_output=True, text=True, cwd=str(tmp_path))
    assert proc.returncode != 0
    assert 'ont' in proc.stderr.lower() or 'pacbio' in proc.stderr.lower()


def test_cli_rejects_both_technologies_at_once(cli_locus, tmp_path):
    proc, _ = run_cli(cli_locus, tmp_path, '--ont', expect_ok=False)
    assert proc.returncode != 0, 'passing --pacbio and --ont should be refused'


def test_cli_version_flag():
    proc = subprocess.run([sys.executable, str(SCRIPT), '--version'],
                          capture_output=True, text=True)
    assert proc.returncode == 0
    assert proc.stdout.strip()


def test_cli_writes_the_blocks_file(cli_locus, tmp_path):
    blocks = tmp_path / 'blocks.bed'
    run_cli(cli_locus, tmp_path, '--output_blocks', blocks)
    assert blocks.exists() and blocks.read_text().strip()


@pytest.fixture(scope='module')
def cli_help():
    proc = subprocess.run([sys.executable, str(SCRIPT), '--help'],
                          capture_output=True, text=True)
    assert proc.returncode == 0
    return proc.stdout


@pytest.mark.parametrize('flag', ['--vcf', '--bam', '--reference', '--chrom',
                                  '--methylation_calls', '--min_mapq',
                                  '--min_base_quality', '--max_allele_length',
                                  '--force'])
def test_documented_flags_are_present(cli_help, flag):
    assert flag in cli_help


@pytest.mark.parametrize('flag', ['--sample', '--llr_thresh', '--error_rate',
                                  '--max_meth_distance'])
def test_cli_exposes_every_constructor_knob(cli_help, flag):
    assert flag in cli_help, f'{flag} is not a command line option'


def test_unknown_flags_are_rejected(cli_locus, tmp_path):
    """Confirms the check above is meaningful: argparse does refuse unknowns."""
    proc, _ = run_cli(cli_locus, tmp_path, '--not_a_real_flag', 'x', expect_ok=False)
    assert proc.returncode != 0
    assert 'unrecognized arguments' in proc.stderr


def test_min_allele_count_default_matches_between_cli_and_constructor():
    """A library caller and a CLI caller must get the same default."""
    import argparse
    import inspect

    from longhap import LongHap, main
    ctor_default = inspect.signature(LongHap.__init__).parameters['min_allele_count'].default
    proc = subprocess.run([sys.executable, str(SCRIPT), '--help'],
                          capture_output=True, text=True)
    assert f'[{ctor_default}]' in proc.stdout, (
        f'help text does not document the constructor default {ctor_default}')


# --------------------------------------------------------------------------- #
# intermediate-file cache
# --------------------------------------------------------------------------- #
CACHE_FLAGS = {
    '--output_transition_matrix': 'transitions',
    '--output_allele_coverage': 'coverage',
    '--output_read_states': 'read_states.json',
    '--output_variant_read_mapping': 'var_reads.json',
    '--output_unphaseable_variants': 'unphaseable',
}


def test_cache_files_are_created(cli_locus, tmp_path):
    extra = []
    for flag, name in CACHE_FLAGS.items():
        extra += [flag, str(tmp_path / name)]
    run_cli(cli_locus, tmp_path, *extra)
    for name in CACHE_FLAGS.values():
        produced = (tmp_path / name).exists() or (tmp_path / (name + '.npz')).exists()
        assert produced, f'{name} was not written in any form'


def test_numpy_cache_path_is_normalised_to_npz(cli_locus, tmp_path, read_phased_vcf):
    """A cache path without a suffix gets ``.npz`` and still round-trips.

    ``np.savez`` always appends ``.npz``, so ``__init__`` normalises the four
    ``--output_*`` numpy paths up front.  The file therefore appears at
    ``<path>.npz``, and -- the part that matters -- the ``os.path.isfile`` reuse
    check in ``infer_variant_transitions`` looks at that same normalised path,
    so a second run hits the cache instead of recomputing.
    """
    target = tmp_path / 'transitions'
    _, first = run_cli(cli_locus, tmp_path, '--output_transition_matrix', target)
    first_records = read_phased_vcf(str(first))

    assert not target.exists(), 'the bare path should not be written'
    assert target.with_suffix('.npz').exists(), (
        'nothing was written at the normalised path either')

    # the reuse check must find the normalised path: same answer, no recompute
    _, second = run_cli(cli_locus, tmp_path, '--output_transition_matrix', target)
    assert read_phased_vcf(str(second)) == first_records


def test_numpy_cache_path_already_ending_in_npz_is_left_alone(cli_locus, tmp_path):
    """An explicit .npz name must not become transitions.npz.npz."""
    target = tmp_path / 'transitions.npz'
    run_cli(cli_locus, tmp_path, '--output_transition_matrix', target)
    assert target.exists()
    assert not (tmp_path / 'transitions.npz.npz').exists()


def test_cache_round_trips_when_paths_already_end_in_npz(cli_locus, tmp_path,
                                                         read_phased_vcf):
    """With explicit .npz names the cache path is exercised end to end."""
    extra = []
    for flag, name in CACHE_FLAGS.items():
        suffix = '.npz' if name in ('transitions', 'coverage', 'unphaseable') else ''
        extra += [flag, str(tmp_path / (name + suffix))]

    _, first = run_cli(cli_locus, tmp_path, *extra)
    first_records = read_phased_vcf(str(first))
    assert (tmp_path / 'transitions.npz').exists()

    # second run must reuse the cache and reproduce the same phasing
    _, second = run_cli(cli_locus, tmp_path, *extra)
    assert read_phased_vcf(str(second)) == first_records


def test_force_recomputes_instead_of_reusing(cli_locus, tmp_path, read_phased_vcf):
    extra = ['--output_transition_matrix', str(tmp_path / 't.npz'),
             '--output_allele_coverage', str(tmp_path / 'c.npz'),
             '--output_read_states', str(tmp_path / 'rs.json'),
             '--output_variant_read_mapping', str(tmp_path / 'vr.json'),
             '--output_unphaseable_variants', str(tmp_path / 'u.npz')]
    _, first = run_cli(cli_locus, tmp_path, *extra)
    baseline = read_phased_vcf(str(first))
    _, second = run_cli(cli_locus, tmp_path, *extra, '--force')
    assert read_phased_vcf(str(second)) == baseline
