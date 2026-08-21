# longhap tests

Hermetic tests for `longhap.py`. Every fixture builds its own reference FASTA,
bgzipped + tabixed VCF and coordinate-sorted BAM in `tmp_path` using pysam, so
the suite needs no reference data, no `example/` files and no external binaries.

## Running

```bash
python -m pytest -q            # from the repo root
python -m pytest tests/test_units.py -q -rxX
```

Requires the runtime dependencies plus pytest:

```bash
pip install -e ".[test]"
```

## Layout

| File | Covers |
|---|---|
| `conftest.py` | Locus builders: reference sequence, variants, read simulator, MM/ML methylation tags, pb-CpG-tools BED |
| `test_fixtures.py` | The builders themselves, checked against pysam |
| `test_units.py` | Pure helpers; the CIGAR walker is checked against `pysam.get_aligned_pairs` |
| `test_end_to_end.py` | Full pipeline on loci with a known truth haplotype |
| `test_semantics.py` | Phasing behaviour: determinism, noise tolerance, switch errors, block breaks, haplotag agreement |
| `test_edge_cases.py` | Degenerate inputs and output plumbing |
| `test_methylation.py` | The methylation-assisted path |
| `test_cli_and_caching.py` | `main()` via subprocess, and the intermediate-file cache |

## Writing a test

`make_locus` takes a callback so variants can be built against the reference it
just generated:

```python
def test_something(make_locus, run_longhap):
    locus = make_locus(lambda seq: alternating_snvs(seq, [400, 600, 800], '010'))
    lh = run_longhap(locus)
    assert len(lh.block_ends) == 1
```

Pass `read_specs=` to control coverage, `genotypes=` to write records longhap
should skip, and `write=False` to stop before the output files are produced.

Two conventions matter:

- **Phasing is defined up to a global flip.** Compare against the truth *and* its
  complement — `agrees_up_to_flip` in `test_semantics.py`.
- **Count switch errors, not mismatches.** `count_switches` measures changes in
  relative phase, which is the metric the benchmarking notebooks use.

## Known-defect tests

Tests marked `pytest.mark.xfail(strict=True)` pin a defect described in
`../longhap_claude_issues.md`; the marker's `reason` names the issue. They fail
today. Each marker sets `strict=True`, so when a fix lands the test starts
passing and that becomes a suite failure — the prompt to remove the marker along
with the memo entry.

Strictness is per-marker rather than a global `xfail_strict` in `pyproject.toml`,
because that file is shared with the `benchmarking/tests` suite.
