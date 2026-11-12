<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->
<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->

# Issue #65: scientific review and maintainer decision

**Prepared for review. No approval is recorded by this document.**

Review target: [PR #69](https://github.com/JoshuaJewell/MetaManifold-WebUI/pull/69).
The numerical implementation is `ec6d5a8ef2630d180b7a8f203ae7dd5bda4d2293`;
list later documentation commits in the actual decision too. Issue #23 and
ordination PR #260 remain excluded.

## Decisions requested

Joshua is asked to assess the [published data-use basis](miseq-data-use.md),
obtain and assess a named independent scientific review, and decide whether
to release a specific revision with a bounded descriptive claim. Data-use
acceptance, scientific approval and release approval are distinct decisions.

The published terms are documented, including their scope inference. A second
permission request is an available response to a specific maintainer concern;
this packet does not claim one was obtained or automatically require one.
A code signature, CI pass, model proof or agent report is not scientific
review. The reviewer should state relevant qualifications, authorship
involvement and conflicts so independence can be assessed.

## What to inspect

| Topic | Source and reviewer task |
| --- | --- |
| Exact summaries | `src/analysis/numeric_policy.jl`, `exact_summaries.jl`, `derived_analyses.jl`: integer/rational preservation, pooling and undefined zero totals. |
| Boundaries | `src/canonical.jl`, `src/server/routes/validated_statistics.jl`: JSON integer syntax, DuckDB streaming, finite limits and immutable persistence. |
| User interpretation | `frontend/src/components/ValidatedStatisticsPanel.tsx`: exact displays, review-pending notice and absence of inference claims. |
| Real-data evidence | Run_A/run_B counts, sample mapping, input/configuration/software hashes and producer results; use pre-taxonomy `seqtab_nochim.csv`. |
| Rights and identity | `miseq-data-use.md`, source manifest and six original FASTQ comparisons. |
| Proof scope | `numeric-contracts.md`: arithmetic model obligations and runtime assumptions; no whole-application or biological correctness proof. |

At the implementation revision, [CI run 42403921664](https://github.com/JoshuaJewell/MetaManifold-WebUI/actions/runs/42403921664)
passed Julia tests (6,337 assertions; five marked broken/skipped), frontend
tests/build, Agda and REUSE. Independent checks passed all 402 real-data
proportions exactly and compared the approximate metrics against Decimal
references. Local and CI matrices have identical hashes. Preserve the CI
evidence before its 7 December 2026 expiry, and recheck CI on the approved head.
That run and the mutation campaign below included a high-precision diversity
layer that is now proposed in a separate pull request.

The six-mutant campaign passed at `419f2fb` after a 499-assertion pristine
baseline. The later minimum-integer rounding correction has eleven new
regressions that failed before and passed after; the final combined policy
and evidence tests passed 161 assertions. Do not claim that the mutation
campaign was repeated on a later revision. Codecov reported 96.86941% patch
coverage and 35 uncovered lines on `ec6d5a8`; neither coverage nor mutation
testing proves that all paths or defects have been covered.

An earlier local differential-route failure also reproduced on untouched
upstream `89be534`, but did not occur in the passing hosted CI. Its cause
remains unexplained and the logs are retained. This contribution changes
no differential implementation or test.

## Independent scientific decision — fill only after review

```text
Reviewer name and relevant qualifications:
Date (UTC):
Project relationship, authorship involvement and conflicts:
Exact commit(s) reviewed:
Evidence archive identifier and SHA-256:

Exact counts / pooling / rational proportions / zero totals: [findings]
JSON / storage / display boundaries and resource limits: [findings]
Sample identity, provenance and biological scope: [findings]
Trusted arithmetic/runtime assumptions and failure modes: [findings]
Local/CI discrepancy, skipped checks and uncovered paths: [assessment]
Reproduction performed (commands, environment, results), or reliance stated:
Required corrections and unresolved questions:

Decision: [APPROVE bounded descriptive scope / CHANGES REQUIRED / DECLINE]
Exact approved claims, limitations and conditions:
Reviewer-authored decision URL or signed record:
```

The proposed claim is exact descriptive totals, richness and rational
proportions on the selected inputs. It excludes hypothesis testing, confidence
intervals, coverage guarantees, taxonomy validation, independent replication
of the three time points, and representative biological performance.

## Joshua's decision — fill only after assessment

```text
Maintainer name:
Date (UTC):
Exact PR head and upstream base reviewed:

Published data-use basis: [ACCEPT / NEED PUBLISHER CLARIFICATION / REJECT]
Reason, attribution location and any specific scope concern:
Named independent scientific decision URL and approved scope:
Remaining conditions and their disposition:

Current-main comparison (behind_by and merge base):
GitHub Verified status for EVERY contribution commit:
Passing checks on the approved head (run URLs):
Durable evidence location and hashes:
CONTRIBUTING requirements reviewed, including fast-forward-only integration:

Release: [APPROVE specified revision and scope / HOLD / DECLINE]
Approved public wording and conditions:
Issue #65: [ELIGIBLE AFTER RECORDED CONDITIONS / KEEP OPEN]
Maintainer-authored decision URL or signed record:
```

Do not populate another person's name, qualification, approval or date.
Scientific approval does not confer publisher rights or release authority;
data-use acceptance does not replace scientific review. Material changes
after a decision require review of their changed scope.

## After real decisions exist

Record their links, retain evidence and recheck upstream main, signatures
and required checks. Follow `CONTRIBUTING.md`; do not rewrite reviewed
history or integrate without the maintainer's decision.

The current validator deliberately requires catalogue status `pending_review`,
and the UI displays that state. Changing only the configuration to `approved`
would fail validation. A future public approval label needs consistent
catalogue, UI, documentation and test changes backed by the actual decisions.
This packet changes no approval state. Close #65 only when its applicable
criteria and recorded decisions are genuinely satisfied.
