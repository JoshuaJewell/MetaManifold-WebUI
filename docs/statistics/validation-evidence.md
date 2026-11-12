<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->
<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->

# Descriptive validation and release gates

Tracks upstream issue #65. The implementation catalogue says `pending_review`.
No fork approval, test pass, proof check or agent-authored document is an
independent statistical review or upstream release approval.

## Reproduce

Use the toolchain in `config/defaults/tool_versions.yml` and the committed
Julia and R lockfiles. Run the normal repository suite, including the server
and pipeline integrations, and `bun run test` / `bun run build` in `frontend`.
The new Julia unit tests are registered in `test/runtests.jl`.

For a generated counts-only feature-by-sample CSV:

```sh
julia --project=. scripts/statistics/biological_validation.jl counts.csv provenance.json result.json
```

The producer reads every CSV field as text and parses counts exactly. Its
arithmetic is checked against committed known answers, not against a second
program at run time. MetaManifold has no Python, so there is no Python oracle:

- Exact totals and reduced proportions: `test/unit/test_exact_summaries.jl`
  compares two tables with the fractions worked by hand in its comments,
  including a count of 2^80 + 17, at zero tolerance.
- Shannon and Gini–Simpson: `test/unit/test_statistics_evidence.jl` checks
  counts 1 and 2 against Gini–Simpson 4/9 exactly and Shannon
  ln 3 − (2/3) ln 2, whose 100 digits come from `bc -l`, with a fixed
  `1e-65` absolute tolerance against 256-bit MPFR. That tolerance is a test
  threshold, not a proved universal error bound.
- High-precision constants: `test/unit/test_high_precision.jl` documents
  where each constant comes from and the `bc` command that reproduces it.

None of these is a run-time check of a real pipeline matrix. To review a real
run, `scripts/statistics/collect_pipeline_evidence.jl` records the producer's
result with the run's provenance: every installed R package against the
lockfile, and the versions and executable hashes of the pipeline tools. CI
does not run it.

Before starting a pipeline run, record its start time with
`export STATISTICS_PIPELINE_STARTED_UNIX=$(date +%s)`. The collector requires
this timestamp and refuses absent, empty or older matrices, so stale output
cannot stand in for a skipped run.

Retain the CSV and its SHA-256, the generated result and check receipts, sample
mapping, original FASTQ hashes, effective configuration, Julia/R/package/tool
versions, source commit and local diff, commands, logs, and lockfile hashes.
Retain this evidence as a CI artefact or review attachment rather than adding
pipeline output, tool archives or build products to the source repository.

## Biological interpretation

The committed MiSeq SOP subset has F3D0, F3D1 and F3D2. The
[source SOP](https://mothur.org/wiki/miseq_sop/) identifies these as time points
from Female 3, not independent biological replicates. It is a small pipeline
exercise, not evidence of performance across representative microbiome studies.
Validate the DADA2 `seqtab_nochim.csv` before taxonomy; `merged.csv` includes
metadata and is not the validator's input format. PR2 is configured in this
bacterial 16S example. This validation makes no taxonomy claim and does not
attempt to establish that reference database's suitability.

Attribution: Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD (2013),
[Development of a dual-index sequencing strategy and curation pipeline for
analyzing amplicon sequence data on the MiSeq Illumina sequencing platform](https://doi.org/10.1128/AEM.01043-13),
Applied and Environmental Microbiology 79(17):5112–5120. Mothur MiSeq SOP
accessed 7 November 2026. The paper links the Mothur-hosted sequences. The
publisher's site-wide notice and upstream's explicit FASTQ annotations give
a published CC BY 4.0 basis. All six decompressed inputs match the original
SOP archive byte-for-byte. The [data-use record](miseq-data-use.md) documents
the sources, attribution and scope inference about the download. This
corrects the earlier incomplete description of the published terms; it does
not supply a new dataset-specific permission or a release approval.

## PR review regressions

The run-route regressions in `test/unit/test_validated_statistics.jl` cover
default and saved table names, typed filters, unknown filter columns, sample
names colliding with projected feature labels, and subgroup names such as
`A` and `AB`. Reserved integer names such as `total` are explicitly rejected
when their count/metadata role is ambiguous. Subgroups use the same underscore-delimited prefix contract as
the existing analysis routes. Tests check the persisted results and selected
sample columns as well as the response. Quoted column identifiers are exercised
through executed text/include/exclude predicates and checked against the exact
selected counts; numeric filters on count columns remain explicitly refused
before the shared floating-point SQL predicate can be used. BIGNUM columns
cannot be silently omitted, and grouped runs persist and retrieve the same
resolved group in their provenance.

CI retains raw Julia coverage even when a test assertion fails, and lists
uncovered descriptive-statistics source lines in the job summary and check
annotations. These reports are diagnostic, not a replacement for a passing
suite; only a passing test run publishes its coverage to Codecov.

`frontend/src/components/AnalysisWorkspace.test.ts` renders the workspace and
descriptive panel, toggles aggregation, and checks the submitted scope,
displayed total, downloaded JSON and rejection of a stale response. Its pinned
React 18 test renderer is a development dependency; the package and lockfile
change supports this interaction test without adding production code.

The simulation now measures absolute and relative errors inside the same
4096-bit precision scope as its reference computation. A regression checks
differences with significant terms more than 256 bits apart under several
ambient precisions, and verifies restoration of MPFR state. New reports record
`comparison_bits`; earlier reports lacking it did not measure the error at the
declared reference precision and must not substitute for the corrected run.

## Required decisions that remain separate

| Gate | Evidence required |
| --- | --- |
| Numerical implementation | Passing known-answer, boundary, adversarial and independent-reference checks at a recorded revision |
| Real pipeline | Retained generated matrix, source/configuration/software provenance and an independent passing receipt |
| Data use | Published CC BY 4.0 basis and attribution documented in `miseq-data-use.md`; record the maintainer's assessment of any specific scope concern |
| Independent scientific review | A named qualified reviewer, scope reviewed, findings, limitations, date and explicit decision |
| Upstream release | JoshuaJewell's explicit approval after the focused PR and required CI pass |

Never populate reviewer or approver fields on someone else's behalf. Keep the
issue open while any applicable gate is outstanding.

Use the [review and decision forms](issue-65-review.md) to record actual
scientific and maintainer decisions for a specific revision.
