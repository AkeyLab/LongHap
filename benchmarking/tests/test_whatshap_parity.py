"""Parity with `whatshap compare`, the only oracle outside this codebase.

Skipped unless whatshap, bcftools and the HG002 chr13 data are all present, so the
rest of the suite stays portable.  This is the check that the variant counting and
the within-block / chromosome-wide switch statistics are not merely self-consistent.
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import SCRIPT

SCRATCH = Path('/scratch/gpfs/AKEY/apfennig/LongHap')
QUERY = (SCRATCH / 'HG002/deepvariant/deepvariant_hifi_hs1.minimap2.hifi.chr13'
                   '.filtered_snps_indels.longhap_meth_phased.vcf.gz')
TRUTH = SCRATCH / 'reference/CHM13v2.0_HG2-T2TQ100-V1.1_stvar.vcf.gz'

WHATSHAP = shutil.which('whatshap') or \
    '/projects/AKEY/akey_vol2/apfennig/micromamba/envs/whatshap/bin/whatshap'
BCFTOOLS = shutil.which('bcftools') or '/projects/AKEY/akey_vol1/software/bin/bcftools'

available = (Path(WHATSHAP).exists() and Path(BCFTOOLS).exists()
             and QUERY.exists() and TRUTH.exists())
pytestmark = pytest.mark.skipif(
    not available, reason='whatshap, bcftools or the HG002 chr13 data is not available')


@pytest.fixture(scope='module')
def normalized(tmp_path_factory):
    """bcftools norm -m -any on both files, as the real invocation does."""
    out = tmp_path_factory.mktemp('norm')
    query, truth = out / 'query.vcf.gz', out / 'truth.vcf.gz'
    subprocess.run([BCFTOOLS, 'norm', '-m', '-any', str(QUERY), '-Oz', '-o', str(query)],
                   check=True, capture_output=True)
    subprocess.run([BCFTOOLS, 'norm', '-r', 'chr13', '-m', '-any', str(TRUTH),
                    '-Oz', '-o', str(truth)], check=True, capture_output=True)
    nops = out / 'query.nops.vcf.gz'
    subprocess.run(f'{BCFTOOLS} annotate -x FORMAT/PS {query} -Oz -o {nops}',
                   shell=True, check=True, capture_output=True)
    return query, truth, nops


def _whatshap(truth, query, out_bed=None, extra=()):
    cmd = [WHATSHAP, 'compare', '--names', 'truth,query', '--ignore-sample-name', *extra]
    if out_bed:
        cmd += ['--switch-error-bed', str(out_bed)]
    cmd += [str(truth), str(query)]
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def _ours(query, truth, out_bed=None, extra=()):
    cmd = [sys.executable, str(SCRIPT), '--vcf', str(query), '--gt_vcf', str(truth), *extra]
    if out_bed:
        cmd += ['--switch_error_bed', str(out_bed)]
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def _stats(text):
    """Map 'label' -> last whitespace-separated token, for lines that have a colon."""
    out = {}
    for line in text.splitlines():
        label, sep, rhs = line.partition(':')
        if sep and rhs.split():
            out[label.strip()] = rhs.split()
    return out


SECTION_DIVIDER = re.compile(r': -{3,}\s*$')


def _section(text, header):
    """One report section, delimited by the next `LABEL: ------` divider.

    Both reports repeat labels across sections -- whatshap prints `switch errors`
    under ALL INTERSECTION BLOCKS and again under LARGEST INTERSECTION BLOCK -- so
    every comparison has to name the section it means.
    """
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if header in line)
    end = len(lines)
    for i in range(start + 1, len(lines)):
        if SECTION_DIVIDER.search(lines[i]):
            end = i
            break
    return '\n'.join(lines[start:end])


def _all_blocks(text):
    """whatshap's pooled statistics, not its largest-block ones."""
    return _stats(_section(text, 'ALL INTERSECTION BLOCKS'))


def test_variant_counts_match(normalized):
    query, truth, nops = normalized
    theirs = _stats(_whatshap(truth, nops))
    ours = _stats(_ours(query, truth))
    for label in ('truth', 'query', 'UNION', 'INTERSECTION'):
        assert ours[label][:3] == theirs[label][:3], label
    assert ours['common heterozygous variants'] == theirs['common heterozygous variants']


def test_chromosome_wide_matches_whatshap_with_ps_stripped(normalized, tmp_path):
    query, truth, nops = normalized
    report = _whatshap(truth, nops)
    theirs = _all_blocks(report)
    ours = _stats(_section(_ours(query, truth), 'CHROMOSOME-WIDE'))
    assert ours['switch errors'] == theirs['switch errors']
    assert ours['phased pairs of variants assessed'] == theirs['phased pairs of variants assessed']
    assert ours['switch/flip decomposition'] == theirs['switch/flip decomposition']
    assert ours['Block-wise Hamming distance'] == theirs['Block-wise Hamming distance']


def test_within_block_matches_whatshap_with_ps_kept(normalized):
    query, truth, _ = normalized
    report = _whatshap(truth, query)
    theirs = _all_blocks(report)
    ours = _stats(_section(_ours(query, truth), 'WITHIN PHASE BLOCKS'))
    assert ours['switch errors'] == theirs['switch errors']
    assert ours['phased pairs of variants assessed'] == theirs['phased pairs of variants assessed']
    assert ours['switch/flip decomposition'] == theirs['switch/flip decomposition']

    full, theirs_pre = _stats(_ours(query, truth)), _stats(report)
    assert (full['non-singleton intersection blocks']
            == theirs_pre['non-singleton intersection blocks'])

    # The per-file block statistics were once computed over the intersection rather
    # than over each file's own phased variants, which whatshap disagreed with.
    for label in ('non-singleton blocks in truth', 'non-singleton blocks in query'):
        assert full[label] == theirs_pre[label], label
    ours_lines, theirs_lines = _ours(query, truth).splitlines(), report.splitlines()
    def covered(lines):
        return [l.split(':')[1].split()[0] for l in lines if '--> covered variants' in l]
    assert covered(ours_lines)[:3] == covered(theirs_lines)[:3]


def test_only_snvs_matches(normalized):
    query, truth, nops = normalized
    theirs = _all_blocks(_whatshap(truth, nops, extra=['--only-snvs']))
    ours = _stats(_section(_ours(query, truth, extra=['--only_snvs']), 'CHROMOSOME-WIDE'))
    assert ours['switch errors'] == theirs['switch errors']
    assert ours['phased pairs of variants assessed'] == theirs['phased pairs of variants assessed']


def test_switch_error_bed_is_identical(normalized, tmp_path):
    """Byte-for-byte on the first three columns, PS honoured on both sides."""
    query, truth, _ = normalized
    theirs_bed, ours_bed = tmp_path / 'theirs.bed', tmp_path / 'ours.bed'
    _whatshap(truth, query, out_bed=theirs_bed)
    _ours(query, truth, out_bed=ours_bed)

    def rows(path):
        return sorted(tuple(line.split('\t')[:3]) for line in
                      path.read_text().splitlines() if line.strip())

    assert rows(theirs_bed) == rows(ours_bed)
    assert len(rows(ours_bed)) > 0


def test_known_chr13_numbers(normalized):
    """Guards the figures this tool has been reporting for HG002 chr13."""
    query, truth, _ = normalized
    ours = _ours(query, truth, extra=['--evaluate_sv'])
    within = _stats(_section(ours, 'WITHIN PHASE BLOCKS'))
    chrwide = _stats(_section(ours, 'CHROMOSOME-WIDE'))
    junction = _stats(_section(ours, 'BETWEEN PHASE BLOCKS'))
    sv = _stats(_section(ours, 'SV PHASING'))

    assert within['switch errors'] == ['111']
    assert within['phased pairs of variants assessed'] == ['96968']
    assert chrwide['switch errors'] == ['218']
    assert junction['junction errors'] == ['107']
    assert junction['block junctions assessed'] == ['191']
    assert sv['non-SNPs in the call set'] == ['17257']
    assert sv['connection errors'] == ['87']
