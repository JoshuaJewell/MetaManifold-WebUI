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
python3 test/references/check_biological_reference.py counts.csv result.json check.json provenance.json
```

The producer reads every CSV field as text and parses counts exactly. The
independent checker rereads the original CSV using a different parser, computes
totals and reduced proportions with Python integers/Fraction, and compares all
values at zero tolerance. It independently computes Shannon and Gini–Simpson
using 100 decimal digits and a fixed `1e-65` absolute tolerance against 256-bit
MPFR. That tolerance is a test threshold, not a proved universal error bound.
The checker refuses Python's optimisation mode, which disables assertions.
It emits a passing receipt only after every comparison passes. When supplied,
the original provenance file is independently hashed and compared too. CI uses
`scripts/statistics/collect_pipeline_evidence.jl` after the integration suite;
this also checks every installed R package against the lockfile and records
the versions and executable hashes of the pipeline tools.

Before starting a pipeline run, record its start time with
`export STATISTICS_PIPELINE_STARTED_UNIX=$(date +%s)`. The collector requires
this timestamp and refuses absent, empty or older matrices. CI records it
before the integration suite so that stale output cannot satisfy a skipped run.

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
accessed 7 November 2026. The paper's data-availability section links to the
Mothur-hosted sequences. Neither that section nor the inspected SOP page
states a separate FASTQ redistribution licence. The repository's code licence does not
by itself establish the FASTQ dataset's redistribution licence. Record the
actual data-use terms or obtain a decision from the responsible data owner
before a redistribution/release claim.

## Required decisions that remain separate

| Gate | Evidence required |
| --- | --- |
| Numerical implementation | Passing known-answer, boundary, adversarial and independent-reference checks at a recorded revision |
| Real pipeline | Retained generated matrix, source/configuration/software provenance and an independent passing receipt |
| Data use | Recorded applicable terms, attribution and any necessary rights approval |
| Independent scientific review | A named qualified reviewer, scope reviewed, findings, limitations, date and explicit decision |
| Upstream release | JoshuaJewell's explicit approval after the focused PR and required CI pass |

Never populate reviewer or approver fields on someone else's behalf. Keep the
issue open while any applicable gate is outstanding.
