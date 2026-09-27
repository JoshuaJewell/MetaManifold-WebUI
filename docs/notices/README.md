<!--
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Run Notices — the user-facing notice model

**Status: design, 2026-09-26.** Nothing here is implemented. It is written
before the code, in the style this repository already uses for
`docs/statistics/method-conditions/*`: where this document and the code
disagree, the document is right and the code is the bug.

**Why it exists.** Today every notice the system produces is a `String` pushed
onto a `Vector{String}`. The strings are individually good — often better than
most software manages — but they carry no category, no severity, no stable
identifier and no statement of whether the numbers may be used. Consequences:

- the UI cannot act on them, so `DangerBanner.tsx` re-derives its own four
  reasons from the config while `AnalysisConfig.jl:1345-1360` derives a
  different list for the backend banner. Two lists, one question, no shared
  source — the drift this repository's own comments warn about when a hook and
  a CI check re-implement the same rule.
- a run whose estimator did not run can be reported safe
  (`docs/audit/2026-09-26-memory-numerics-warning-audit.md` §W1), because
  nothing in the notice says "this blocks export".
- a check that is a stub and a check that found nothing look identical
  (`has_batch_confounding => false`, §N3).
- the log and the result cannot be correlated: the log has a timestamp and the
  manifest has a list, and there is no key between them.

**Non-goals.** This is not an error-handling framework, not a replacement for
exceptions, and not a logging library. Exceptions still mean "we stopped";
`@warn`/`@error` still mean "operator, look at the log". A `RunNotice` means
*"someone reading these results has to know this"*, and it travels **with the
result**, not only with the log.

---

## 1. The six categories

The category is decided by *the question the notice answers*, not by where it
is raised. Getting this wrong is how "we did not check" gets read as "we
checked and found nothing".

| Category | Question it answers | Typical examples |
|---|---|---|
| `technical_resource` | Did the machine have the room, the time, the precision? | allocation churn, peak memory, a benchmark threshold, a check that is **not implemented**, R warnings raised outside the handler |
| `dependency_environment` | Was the toolchain there and reachable? | R absent, `vegan` absent, `MASS` absent, R runtime busy, CI produced no verdict |
| `data_quality` | Is the input table what we think it is? | zero variance, library-size outliers, low prevalence/abundance, all-zero samples or taxa, samples dropped below depth |
| `method_limitation` | Does the chosen method mean what a reader will think it means? | data-derived dendrogram basis, non-isometric balance weights, `rarefy`, an inert configuration parameter, `theta` at a boundary |
| `result_safety` | May these numbers be read as a finding? | estimation `:not_run`, NaN-filled ordination coordinates, healing that rewrote values, SBP p-hacking guard |
| `fatal` | Did we stop? | hard-stop with a DANGER banner, refusal of a refused dispersion method, refusal of an unrooted tree |

Two rules that fall out of the table and are worth stating as invariants:

- **A diagnostic that did not run is `technical_resource`, never
  `data_quality`.** Its summary must say it did not run.
- **"Nothing was produced" is `result_safety` (or `fatal`), never silence.** A
  run with no statistics is not a safe run.

## 2. Severity

| Severity | Meaning | UI treatment |
|---|---|---|
| `info` | recorded for provenance; nobody needs to act | collapsed, one line |
| `notice` | worth knowing; changes how you read a number | collapsed, one line |
| `warning` | interpretation required before use | expanded, amber |
| `danger` | do not present these numbers as a finding without saying this | expanded, red, banner on every figure |
| `fatal` | nothing was produced | replaces the result view |

**Invariants** (these are the assertions in `test/unit/test_notices.jl`):

1. `result_usable == false` ⟹ `blocks_export == true`.
2. `severity == fatal` ⟹ `result_usable == false` **and** the result object
   holds no statistics.
3. `category == fatal` ⟺ `severity == fatal`.
4. `category == technical_resource` ⟹ `result_usable == true`. (The machine
   struggling is never a reason to distrust the numbers.)
5. `category == result_safety` ⟹ `blocks_export == true`.
6. Every notice has a non-empty `summary`, `detail` and `id`; `action` may be
   empty only when `severity ∈ {info, notice}`.
7. Every `id` emitted anywhere in the codebase appears in
   `docs/notices/catalogue.md`, and every `id` in the catalogue is emitted
   somewhere. This is the same drift guard `PipelineLog._LOG_REGISTRY` already
   enforces for log files (`src/core/log.jl:77-84`, `test/unit/test_log.jl`),
   applied to notices.

## 3. The field model

```julia
module Notices

using Dates, OrderedCollections

export RunNotice, NoticeCategory, NoticeSeverity,
       TECHNICAL_RESOURCE, DEPENDENCY_ENVIRONMENT, DATA_QUALITY,
       METHOD_LIMITATION, RESULT_SAFETY, FATAL_CATEGORY,
       emit!, notices, blocks_export?, result_usable?

@enum NoticeCategory technical_resource dependency_environment data_quality
                     method_limitation result_safety fatal_category

@enum NoticeSeverity info notice warning danger fatal

"""
    RunNotice

One thing a reader of these results has to know. See docs/notices/README.md.
`rev` is the catalogue revision of this id's wording, so a stored notice can be
rendered against the catalogue that was current when it was emitted.
"""
struct RunNotice
    id             :: String                  # stable, dotted, see §4
    category       :: NoticeCategory
    severity       :: NoticeSeverity
    summary        :: String                  # ≤ 120 chars, one line, UI
    detail         :: String                  # paragraphs, no length limit
    action         :: String                  # what to do; "" if nothing
    result_usable  :: Bool
    blocks_export  :: Bool
    log_ref        :: String                  # see §5
    data           :: OrderedDict{String,Any} # structured context, JSON-safe
    rev            :: Int
end
```

The nine fields the brief asks for are `category`, `severity`, `id`, `summary`,
`detail`, `action`, `result_usable`, `blocks_export` and `log_ref`. The two
additions:

- **`data`** — structured context (`bytes`, `taxa`, `threshold`, `count`,
  `package`, `waited_seconds`). Without it every consumer has to parse the
  `summary` to find the number, and the summary has to be rewritten to change
  a threshold. Every notice carries its numbers here, and the summary
  interpolates them.
- **`rev`** — the catalogue revision. A manifest written last year has to
  render against last year's wording, or a reworded notice silently changes
  the meaning of a stored result.

### Serialisation

One shape, three targets (manifest, HTTP response, DOI bundle). JSON:

```json
{
  "id": "ilr.dendrogram.data_derived",
  "rev": 1,
  "category": "method_limitation",
  "severity": "warning",
  "summary": "The balance dendrogram was built from the same table that is tested.",
  "detail": "The variation matrix Var(log(x_i/x_j)) ... no sample metadata, so it cannot see group labels, but it is not independent of the data.",
  "action": "Pre-register the dendrogram method, or use a phylogenetic or SBP basis chosen without reference to this table.",
  "result_usable": true,
  "blocks_export": false,
  "banner_required": true,
  "log_ref": "run=7c1d… seq=17 origin=src/analysis/ilr_basis.jl:1265",
  "data": { "taxa": 1204, "method": "ward", "variation_matrix_bytes": 5798304 }
}
```

`banner_required` is **derived**, not stored: `severity ∈ {danger, fatal} ∨
(category == method_limitation ∧ severity == warning)`. The repository already
promises "Bannered in every figure and DOI bundle" for dangerous
configurations; deriving the flag means the UI and the backend cannot disagree
about which notices get a banner — the failure mode described at the top of
this document.

The Nickel contract belongs beside the other schemas in `config/schemas/`,
following `analysis_config.ncl` and `doi_publication.ncl`.

## 4. Stable identifiers

Rules, so that the set is checkable and grep-able:

```
<area>.<subject>.<condition>
```

- `<area>` is the module that owns the semantics:
  `config`, `exec`, `estimation`, `ilr`, `scaling`, `r`, `analysis`, `bench`,
  `provenance`, `pipeline`, `ci`.
- `<subject>.<condition>` is lowercase `snake_case`, no digits that encode
  values, no free text, no interpolation. `ilr.sbp.p_hacking_guard`, not
  `ilr.sbp.3_matrices_tried`.
- An `id` is **immutable**. To change what a notice means, add a new `id` and
  retire the old one (retired ids stay in the catalogue with
  `status: retired` so old manifests still render).
- Values never live in the `id`; they live in `data`. `exec.heal.nan_inf`
  carries `{"nan_count": 3, "inf_count": 0, "epsilon": 1e-6}` — so the
  catalogue has one entry, not one per count.

## 5. The technical log reference

`log_ref` answers *"where do I look for the details?"* — the question a
scientist asks after reading the summary. The repository already has the
precedent: `PipelineLog.CMD_MARKER` exists so that
`grep '^\[MetaManifold\] cmd: '` recovers every command a run issued
(`src/core/log.jl:53-59`).

Format:

```
run=<run_id> seq=<n> origin=<file>:<line> [logger=<key>] [log=<relative path>]
```

- `run_id` — a UUID per `prepare_analysis_table`/`run_analysis` pair *and* per
  HTTP request that produces a chart. **New:** the server has no correlation
  id today.
- `seq` — a per-run counter, so two healings in one run are distinguishable
  and the log order is recoverable.
- `origin` — the emission site. Recorded by the emitter, not derived, because
  `@__FILE__`/`@__LINE__` at the call site is the only thing that survives a
  refactor of the message text.
- `log` — for pipeline stages, the registry-relative path from
  `PipelineLog._TOOL_LOG_FILES`, so `grep "$run_id" combined_pipeline.log`
  works.

The convention: `grep "$run_id"` in the process log returns every log line the
run produced, in order. That is one `@info`-with-`run_id` change in the emitter
and it makes every notice traceable.

## 6. Emission

Notices are emitted **alongside** the existing free-text warnings, not instead
of them, until every consumer has migrated. Concretely:

```julia
push!(warnings, "Low prevalence/abundance: $(n_low) taxa below min_prevalence=$(mp) — will be filtered")
emit!(sink, RunNotice(
    id = "exec.filter.low_prevalence", rev = 1,
    category = DATA_QUALITY, severity = :warning,
    summary = "$(n_low) of $(n_taxa) taxa are below min_prevalence=$(mp) and were filtered.",
    detail  = "…",
    action  = "…",
    result_usable = true, blocks_export = false,
    log_ref = log_ref(run, origin),
    data    = OrderedDict("n_filtered" => n_low, "threshold" => mp)))
```

Two rules for the migration:

- **No notice is emitted without its free-text twin, and vice versa**, until
  the drift test in §2(7) has been green for one release.
- **Nothing is removed from the manifest.** `diagnostics.warnings` stays;
  `diagnostics.notices` is added. A manifest that already exists must still
  render.

Emission points, in the order they are worth doing:

| where | what becomes notices |
|---|---|
| `Execution.run_analysis` | `estimation.not_run`, and the `is_dangerous` fix (§W1) |
| `Execution.self_diagnostics` | all eight `warnings` strings, plus `exec.batch_confounding.not_implemented` |
| `Execution.safe_self_healing` / the healing block | `exec.heal.*`, with positions (§N2) |
| `ILRBasis.ilr_transform` | the four `warnings` and the `dangers` |
| `Scaling.*_factors` | the five scaling warnings |
| `Analysis._ensure_r`, `run_nmds`, `run_permanova` | `r.vegan_unavailable`, `r.runtime_busy`, `ord.*.not_run` |
| `Estimation._not_run` / `_assemble` | `estimation.not_run`, `estimation.feature_failed` |
| HTTP routes | attach the response's notice list to the JSON body, so a degraded chart says so in the payload and not only in the log |
| `bench/ilr_bases/*` | `bench.allocation_churn`, `bench.peak_rss_growth` with the phase in `data` |

## 7. What the UI does with them

The contract is small, and it is the whole point of the model:

1. Group by `category`, order by `severity` descending.
2. Render `summary` in the list; `detail` and `action` on expand; `log_ref` as
   a copyable string.
3. **If any notice has `blocks_export`, disable the export/DOI action and say
   which notice blocked it.** Today nothing can ask that question.
4. **If any notice has `result_usable == false`, the numbers are not rendered
   as a results table.** They are rendered as a statement about the run.
5. Every notice with `banner_required` is rendered on every figure, as
   `DangerBanner` already promises for dangerous configurations.
6. `DangerBanner.tsx` stops deriving its own list from the config and renders
   the `config.*` notices the backend emitted — one source, both surfaces.

## 8. What this model deliberately does not do

- **It does not decide policy at emission time.** The category and severity are
  constants in the catalogue; the emitter supplies the data. Otherwise the same
  condition gets classified differently in two places, which is the drift this
  is meant to remove.
- **It does not suppress anything.** The brief's "do not suppress warnings
  globally" is satisfied by construction: notices are additive, the free text
  stays, and `severity = info` means "recorded", not "hidden". The audit found
  no global suppression today and this model adds none.
- **It does not replace the DANGER banner.** The banner is the human-readable
  rendering of the `danger` notices plus the config's own acknowledgements. One
  list, rendered two ways.
- **It does not carry results.** A notice describes a run, not a feature.
  Per-feature statuses stay where they are (`status = "failed"` with a `note`
  per row in `Estimation`), and `estimation.feature_failed` is the *aggregate*
  notice that points at them.

## 9. Open questions for the owner

1. **Does `result_usable == false` block export, or require an
   acknowledgement?** The model currently blocks. The repository already has
   `DANGER_ACK_TOKEN` for *configuration* risk; extending the same mechanism to
   *result* risk is a decision, not a default.
2. **Are `technical_resource` notices in the manifest, or only in the log?**
   They are provenance-relevant (a run that peaked at 4 GiB is a fact about the
   run) but they are not science. The catalogue marks each with
   `manifest: yes/no`.
3. **Who owns the `r.*` notices' wording** — MetaManifold, or the standalone
   runtime package that would eventually own the R hand-off? Recommendation:
   MetaManifold, because "vegan is not installed" is a sentence about this
   product's results.

---

## Files

| File | Contents |
|---|---|
| `docs/notices/README.md` | this document — the model |
| `docs/notices/catalogue.md` | every notice id, with category, severity, wording and the audit finding it came from |
