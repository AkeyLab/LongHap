"""Fixtures for the evaluate_phasing tests.

cyvcf2 reads plain uncompressed VCFs, so every fixture is a text file written to
``tmp_path``.  No bgzip, no tabix, no external tools -- the hermetic part of the
suite needs nothing but pytest and cyvcf2.
"""
import gzip
import subprocess
import sys
from pathlib import Path

import pytest

MODULE_DIR = Path(__file__).resolve().parent.parent
SCRIPT = MODULE_DIR / 'evaluate_phasing.py'
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))


VCF_HEADER = """##fileformat=VCFv4.2
{contigs}##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}
"""


def _record(chrom, pos, ref, alt, calls, af=None):
    """One VCF line. ``calls`` is a list of (gt, ps) per sample, ps may be None.

    ``af`` of None writes no INFO/AF, which is how a variant absent from the
    reference panel arrives.
    """
    fmt = 'GT:PS' if any(ps is not None for _, ps in calls) else 'GT'
    fields = []
    for gt, ps in calls:
        fields.append(f'{gt}:{ps}' if 'PS' in fmt else gt)
        if 'PS' in fmt and ps is None:
            fields[-1] = f'{gt}:.'
    info = '.' if af is None else 'AF=' + (','.join(str(a) for a in af)
                                           if isinstance(af, (tuple, list)) else str(af))
    return '\t'.join([chrom, str(pos), '.', ref, alt, '.', 'PASS', info, fmt] + fields)


@pytest.fixture
def write_vcf(tmp_path):
    """Build a VCF from ``(chrom, pos, ref, alt, gt, ps)`` tuples.

    ``gt`` is written verbatim ('0|1', '1/0', '.|1', '1|1'), so a test states the
    genotype exactly as the file would carry it.  ``ps`` of None omits the tag.
    For several samples pass ``gt``/``ps`` as tuples of per-sample values.

    A record may carry a seventh element, the INFO/AF value; omitting it writes no
    AF at all, which is how a variant missing from the reference panel arrives.
    """
    def _write(name, records, samples=('S1',), contigs=('chr1',)):
        path = tmp_path / name
        contig_lines = ''.join(f'##contig=<ID={c},length=100000000>\n' for c in contigs)
        lines = [VCF_HEADER.format(contigs=contig_lines, samples='\t'.join(samples))]
        for record in records:
            chrom, pos, ref, alt, gt, ps = record[:6]
            af = record[6] if len(record) > 6 else None
            gts = gt if isinstance(gt, tuple) else (gt,) * len(samples)
            pss = ps if isinstance(ps, tuple) else (ps,) * len(samples)
            lines.append(_record(chrom, pos, ref, alt, list(zip(gts, pss)), af) + '\n')
        path.write_text(''.join(lines))
        return str(path)
    return _write


@pytest.fixture
def write_bed(tmp_path):
    """Build a 4-column exon BED; ``name`` ending in .gz is gzipped."""
    def _write(name, rows):
        path = tmp_path / name
        text = ''.join('\t'.join(str(f) for f in row) + '\n' for row in rows)
        if name.endswith('.gz'):
            path.write_bytes(gzip.compress(text.encode()))
        else:
            path.write_text(text)
        return str(path)
    return _write


@pytest.fixture
def run_cli(tmp_path):
    """Invoke the script as a subprocess; returns the CompletedProcess."""
    def _run(*args, expect_ok=True):
        proc = subprocess.run([sys.executable, str(SCRIPT), *map(str, args)],
                              capture_output=True, text=True)
        if expect_ok:
            assert proc.returncode == 0, f'exit {proc.returncode}\n{proc.stdout}\n{proc.stderr}'
        return proc
    return _run
