<!--
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Audit — memory, numerical, backend and warning boundaries (2026-09-26)

**Audit type:** read-and-trace audit of `src/analysis/*`, `src/core/r_runtime.jl`,
`src/core/provenance.jl`, `src/analysis/analysis.jl`, `bench/ilr_bases/*` and the CI
state, against the ten questions below.
**Auditor:** Agentic audit. Every finding cites `file:line`. Arithmetic is in
Appendix A so each number can be re-derived.
**Revision audited:** `main` @ `4a848da` (branch `arena/01a0de46-metamanifold-webui`,
clean tree at audit time).

## What I could and could not run

| | |
|---|---|
| Static read of all of `src/analysis/**`, `bench/ilr_bases/**`, `renv.lock` | done |
| `gh run list` / `gh api …/jobs` (CI state) | done — see **B1** |
| `julia --project=. -e 'using Pkg; Pkg.test()'` | **not run** — no Julia toolchain in the audit sandbox, and `julialang.org` is unreachable from it |
| `R`/`vegan`/`MASS` availability | **not verified** — no R in the sandbox |

So every memory number below is *derived from the source*, not measured. Each is
labelled **DERIVED**. The one thing that would turn them into measurements is
`bench/ilr_bases/benchmark.jl`, which already reports `allocated_bytes` and
peak-RSS growth per workload; it has never run in CI (B1), and `bench/results/`
is git-ignored, so **no recorded ILR benchmark output exists anywhere in the
repository.**

## Verdict

The statistical layer is unusually careful: refusal rather than substitution is
enforced in `Estimation`, `Scaling`, `ILRBasis` and `AnalysisConfig`, and the
conditions documents are written before the code. The audit found **no case
where a refused method is silently replaced by another**.

What it did find falls into three groups:

1. **Two ordering defects in the result-safety layer** (W1, B3) that let a run
   with *no statistics at all* be reported as safe, and let one R-dependent
   analysis hand back NaN-filled coordinates while a sibling correctly reports
   "not run". Both are cheap to fix and both are the kind of thing this
   repository already refuses everywhere else.
2. **One real algorithmic regression the new engine has already solved elsewhere**
   (M1): the *default* Helmert ILR basis is O(D²·n) in allocation, while the
   issue-#20 engine that sits beside it is O(D·n) and is already proven equal to
   it at 1e-12 by an existing test.
3. **A set of avoidable duplicate allocations** (M2–M7) that are individually
   modest and collectively the difference between "10 000 taxa is fine" and
   "10 000 taxa is 8 GB of churn".

**Blocker, external to the code:** CI has been `startup_failure` with zero jobs
since 2026-09-25 22:19 UTC, including on the merge commit that landed the ILR
bases (#75). There is **no recorded test verdict on the code being audited**.
The preconditions the owner set — *current tests pass, current CI is
operational, reference outputs captured, parity criteria defined* — are not met,
so **nothing in this audit recommends a new numerical backend.**

---

# Part 0 — The warning taxonomy

Every notice the system can emit is classified into exactly one of six
categories. The categories are chosen so that the *action* follows from the
category, and so that "the check did not run" is never in the same box as "the
check ran and found nothing".

| Category | Question it answers | Default severity | Result usable? | Blocks export? |
|---|---|---|---|---|
| `technical_resource` | Did the machine have the room/time? | `info`/`notice` | yes | no |
| `dependency_environment` | Was the toolchain there? | `warning`/`error` | only if something ran | yes if nothing ran |
| `data_quality` | Is the input table what we think it is? | `warning` | yes, with interpretation | no |
| `method_limitation` | Does the chosen method mean what the reader will think it means? | `warning` | yes, with interpretation | no (banner instead) |
| `result_safety` | May these numbers be read as a finding? | `danger` | **no, not as-is** | yes |
| `fatal` | Did we stop? | `fatal` | no — nothing was produced | yes |

Two rules fall out of the table and both are violated today:

- **A check that did not run is `technical_resource`, not `data_quality`.**
  `check_batch_confounding` (Execution.jl:500-518) always returns
  `has_batch_confounding => false`; that is `technical_resource` ("this check is
  a stub") and must never render as "no confounding detected" (**N3**).
- **A result with no statistics is `result_safety`, and it blocks export even
  when nothing "went wrong".** A run whose estimator returned `:not_run` is not
  a safe run (**W1**).

---

# Part 1 — Findings

Severity is the notice severity the finding *should* carry. IDs are stable
within this audit (`M` memory, `N` numerical, `B` backend, `W` warning
boundary).

## Memory

### M1 — The default Helmert ILR basis is O(D²·n) in allocation, and the O(D·n) replacement is already in the tree `severity: notice · category: technical_resource`

`Execution.jl:1188` — `mean_first_i = mean(log_col[1:i])` inside
`for i in 1:(n_taxa-1)`, inside `for j in 1:size(clr_table, 2)`.

`log_col[1:i]` is a fresh `Vector{Float64}` of length `i`, allocated and summed
once per balance per sample. Σᵢ 8i ≈ 4·D² bytes per sample, so the whole
transform allocates **4·n·D·(D−1) bytes**:

| taxa × samples | allocation per run (DERIVED) |
|---|---|
| 1 000 × 20 | 80 MB |
| 1 500 × 40 (the regression-gate workload) | 360 MB |
| 10 000 × 20 (the benchmark's largest size) | **8.0 GB** |

The comment at `Execution.jl:1175-1177` says the default loop is "left exactly
as it was so that no default-basis result moves", which was the right call for
the PR that changed everything around it. But the engine it is protecting
against drift from is in the same repository, is pure Julia, and is already
asserted equal:

- `ILRBasis.comb_tree` is documented as *"the comb (caterpillar) tree whose
  balances are the Helmert default"* (`ilr_basis.jl:153-160`);
- `test/unit/test_ilr_basis.jl:195-208` (`comb-is-helmert`) computes the default
  basis through `prepare_analysis_table` and through `comb_tree`, and asserts
  agreement to **1e-12**;
- `ILRBasis.tree_balances` is one post-order pass per sample — O(D) time and
  O(D) memory per sample, no dense basis (`ilr_basis.jl:1078-1105`).

So the swap has an existing parity test and an existing tolerance. It is gated
in Part 5, not proposed for today.

**Change (when the gate opens):** route `ilr_basis = "default"` through
`ILRBasis.tree_balances(ILRBasis.comb_tree(taxa), counts_after_zero, ones(D), ones(D-1))`,
keeping the `any(col .<= 0)` refusal. **Belongs in:** MetaManifold.

### M2 — The balance dendrogram peaks at ~2× the size it reports, so both thresholds are off by 2× `severity: warning · category: technical_resource`

`ilr_basis.jl:1250` computes `need = 8·(D·(D−1)÷2)` and compares it with
`DENDROGRAM_WARN_BYTES` (1 GiB) and `DENDROGRAM_REFUSE_BYTES` (2 GiB). That is
the size of **one** condensed variation vector. Two are live at once:

- `variation_condensed` returns `tau` (`ilr_basis.jl:773`);
- `hclust_r` begins `diss = Float64.(diss_in)` (`ilr_basis.jl:803`). Broadcast
  always allocates a new array, so `diss` is a second full copy of a
  `Vector{Float64}` that is still referenced by the caller's `tau`.

| taxa | reported `variation_matrix_bytes` | actual peak (DERIVED) | ratio |
|---|---|---|---|
| 10 000 | 381 MiB | 764 MiB | 2.00× |
| 16 385 (first D that warns) | 1.00 GiB | 2.00 GiB | 2.00× |
| 23 171 (first D that is refused) | 2.00 GiB | **4.00 GiB** | 2.00× |

The refusal message (`ilr_basis.jl:1251`) says *"above the 2 GiB limit"* about a
run that would have peaked at 4 GiB.

**Change (low risk, two lines):** `hclust_r` may mutate `diss_in` in place —
nothing else reads `tau` afterwards — so either take ownership
(`diss = diss_in isa Vector{Float64} ? diss_in : Float64.(diss_in)`) or document
the copy and budget `2·need` in both thresholds. The second option alone is a
one-constant change and is the safer first step.
**Belongs in:** MetaManifold.

### M3 — The ILR branch builds a D×n CLR table and never reads a single value from it `severity: notice · category: technical_resource`

`Execution.jl:1154-1180`:

```
clr_table = similar(counts_after_zero)   # 1154  — D×n Float64
…  clr_table[:, j] = log_col .- mean(log_col)   # 1161  — written
n_taxa = size(clr_table, 1)              # 1164  — read: dimensions only
ilr_table = zeros(n_taxa - 1, size(clr_table, 2))   # 1179 — dimensions only
for j in 1:size(clr_table, 2)            # 1180  — dimensions only
```

`clr_table` is read for its **shape** and never for its contents. Both the
default branch (which recomputes `log.(counts_after_zero[:, j])` at line 1184)
and the engine branch (which takes `counts_after_zero`) ignore it. It costs one
`D×n` allocation plus an `O(D·n)` log-and-centre pass.

The loop is not pure waste: its `any(col .<= 0)` check at line 1158 is the
refusal that guarantees positivity, and that must survive.

**Change (low risk):** replace the loop with a validation-only pass over the
columns (`log_col = log.(@view counts_after_zero[:, j])`; throw on any
non-positive; discard) and take `n_taxa` from `counts_after_zero`. Identical
refusal behaviour, one fewer `D×n` matrix.
**Belongs in:** MetaManifold.

### M4 — `prepare_analysis_table` holds three to four copies of the largest object `severity: notice · category: technical_resource`

| line | object | needed? |
|---|---|---|
| 840 | `filtered_counts = copy(counts)` | **no** — only read to compute prevalence/abundance, then replaced at 858 |
| 858 | `filtered_counts = filtered_counts[keep_taxa, :]` | yes |
| 868 | `filtered_counts = filtered_counts[top_indices, :]` | yes, but composable with 858 into one slice |
| 1059 | `prepared = copy(counts_after_zero)` | only for `clr` (writes in place at 1147) and `rarefy` (1233); every other branch reassigns `prepared` before reading it |
| 1161 | `clr_table` (see M3) | no |

For a 20 000-taxon input filtered to 10 000 taxa over 200 samples the peak drops
from ~107 MiB to ~76 MiB by removing the dead copy, the composable slice and
`clr_table` (DERIVED, Appendix A).

**Change (low risk):** drop the `copy` at 840 (nothing mutates `counts` before
the slice); compose the two row filters into one index vector before slicing
once; move the `copy` at 1059 into the `clr` and `rarefy` branches only.
**Belongs in:** MetaManifold.

### M5 — The SBP path materialises its source three times, and its "streaming" reader allocates per cell `severity: warning · category: technical_resource`

`ilr_basis.jl:1228-1231`:

```
bytes = _read_source(sbp_path, …)                    # the whole file
source_sha = bytes2hex(sha256(bytes))
(sbp_taxa, ids, W) = parse_sbp(String(copy(bytes)))  # a second whole-file copy
```

`String(copy(bytes))` doubles the file in memory. The `copy` is presumably
defensive against `String(::Vector{UInt8})` taking ownership of the buffer —
but `bytes` is not read after the `sha256`, so `String(bytes)` is safe here and
halves the peak. `W` is then `Matrix{Int8}` of `D × (D−1)`.

The reader itself keeps only one row at a time (`ilr_basis.jl:474-479`), which
the docstring is right to advertise — but retention is not allocation. Pass 2
(`ilr_basis.jl:545-570`) allocates a `String` per cell (`String(take!(field))`)
and then a `strip` result per cell, and `parse_sbp` runs the whole reader twice
over the text (pass 1 for shape, pass 2 for values).

| taxa | file text | `String(copy(…))` | `W` (Int8) | peak (DERIVED) | cells |
|---|---|---|---|---|---|
| 1 000 | 1.9 MiB | 1.9 MiB | 1.0 MiB | 5 MiB | 1.0e6 |
| 10 000 | 191 MiB | 191 MiB | 95 MiB | **477 MiB** | 1.0e8 |

At 10 000 taxa that is ~10⁸ cell conversions on top of the 477 MiB, and the file
format itself is 191 MiB on disk before any of it. The `Matrix{Int8}` is the
right representation; the two text copies and the per-cell `String`s are not.

**Change (measure first, then fix):** (a) `String(bytes)` instead of
`String(copy(bytes))` at 1221 and 1231 — one line each, halves the peak;
(b) hash the file in chunks while streaming it, so no `bytes` array exists at
all; (c) parse the CSV from an `IO` into `W` without building a `String` per
cell. (b) and (c) are real work with their own test matrices.
**Belongs in:** (a) MetaManifold; (b) and (c) a standalone runtime package
(streaming readers with a bounded-allocation contract), or upstream if one
exists.

### M6 — DuckDB already is the on-disk intermediate representation; the Julia hand-off duplicates its result `severity: notice · category: technical_resource`

`analysis.jl:108-124` (`filtered_counts`):

```
result = DataFrame(DBInterface.execute(con, sql, where_params))   # F×S, fully materialised
mat    = zeros(Float64, n_samples, n_features)                   # S×F again
for (j, row) in enumerate(eachrow(result))
    for (i, col) in enumerate(sample_cols)
        mat[i, j] = ismissing(row[Symbol(col)]) ? 0.0 : Float64(row[Symbol(col)])
```

DuckDB has already done the filtering on disk. The result is then held **twice**
(a `DataFrame` and a transposed `Matrix`) and filled through `O(F·S)`
`DataFrameRow` name lookups — a hash per cell, 5 000 000 of them at
10 000 × 500.

This answers question 5 directly: the peak-memory problem in the analysis path
is **not** that the data lives in RAM instead of on disk. It is that the
on-disk engine's output is copied and then re-indexed cell by cell.

**Change (low risk):** take the result column-wise —
`for (i, col) in enumerate(sample_cols); mat[i, :] = coalesce.(Tables.getcolumn(result, Symbol(col)), 0.0); end`
— one name lookup per *column* instead of per cell, and a typed vector instead
of a boxed `Any`. Same numbers.
**Belongs in:** MetaManifold.

### M7 — `part_weights` computes all three weight kinds to return one `severity: info · category: technical_resource`

`ilr_basis.jl:976-1000`: `gm`, `an` and `en` are all computed (three passes over
`X`, three `D`-vectors, four `n`-vectors of temporaries per row) and then one of
five combinations is selected. For `part_weights_kind = "gm_counts"` — the
default alternative to `uniform` — two thirds of that work is discarded.

**Change (low risk):** branch to compute only the selected kind. **Belongs in:**
MetaManifold. Do it after M1–M4, not before; it is cosmetic by comparison.

### M8 — `variation_condensed` walks its cache array against the grain `severity: info · category: technical_resource`

`ilr_basis.jl:761` allocates `C = Matrix{Float64}(undef, n, D)` and then fills
it as `C[s, i]` with `i` outer and `s` inner — i.e. writing down a column of a
column-major array, one cache line per element. The later `tau` loop reads
`C[s, i]` and `C[s, j]` over `s`, two strided streams.

Transposing to `Matrix{Float64}(undef, D, n)` with `C[i, s]` makes both loops
contiguous. Footprint is unchanged; this is a locality fix only, worth doing
only if the profiling from Part 5 says the kernel matters.
**Belongs in:** MetaManifold.

## Numerical

### N1 — NaN/Inf healing writes values on the wrong scale `severity: danger · category: result_safety`

`Execution.jl:600-616`:

```
isnan(healed[i])  →  healed[i] = epsilon          # e.g. 1e-6
isinf(healed[i])  →  healed[i] = log(1/epsilon)   # e.g. +13.8   (positive Inf)
                     healed[i] = log(epsilon)     # e.g. −13.8   (negative Inf)
```

The `Inf` branch treats `prepared` as a log scale; the `NaN` branch treats it as
a linear scale. Both branches run on the same table. Consequences:

- On `relative`/`presence_absence`/`tss`-as-proportions, whose values lie in
  [0, 1], a `+Inf` becomes **13.8** — a number the transform cannot produce.
- On `clr`/`ilr`, a `NaN` becomes **1e-6**, which is a near-zero *log-ratio*
  rather than the near-zero *linear* value the same epsilon would mean on a
  proportions table.

Two arbitrary sentinel values, chosen by an inconsistency, written into the
table that is then fitted. The counts are recorded (`Execution.jl:1306`) but not the values, so nothing downstream can tell.

**Change (low risk, and a decision the owner should make explicitly):** pick the
rule per transform scale and record it. The honest options are (a) refuse
(`drop_policy = REFUSE` already exists for exactly this class of situation) or
(b) replace with the transform's own agreed sentinel and record, per healing,
which positions changed. Option (b) needs a mask or a position list; see N2.
**Belongs in:** MetaManifold — this is a science decision, not a runtime one.

### N2 — Self-healing mutates data and records only a count `severity: warning · category: result_safety`

Question 7, answered: **self-healing changes data, not metadata.**

| healing | what changes |
|---|---|
| `heal_nan_inf` (`Execution.jl:599`) | rewrites values **in `prepared`**, the table that is fitted |
| `heal_all_zero_samples` / `heal_all_zero_taxa` with `IMPUTE` (`Execution.jl:654`, `675`) | writes `epsilon` into `filtered_counts`; for all-zero *samples* this happens **before** the transform, so it changes every transform's input |
| … with `DROP` | removes rows/columns, so the result table has a different shape **and** the surviving taxa/samples are renumbered |

All three are recorded — `diagnostics.healings`, and again in
`provenance["diagnostics"]["healings"]` (`Execution.jl:1676-1679`). What is
recorded is a **sentence with a count**:
`"Healed 3 NaN and 0 Inf with epsilon=1e-6"`. There is no record of *which*
cells, so a healed value is indistinguishable from a measured one in the result
table, and two runs with the same input and the same count can differ.

**Change (low risk, additive):** alongside each healing entry, record a compact
fingerprint of what changed — the affected `(row, col)` positions when they are
few, otherwise a count plus a SHA-256 of the position bitmap. The prepared-table
hash (`Execution.jl:740`) already exists as a precedent for fingerprinting a
matrix. **Belongs in:** MetaManifold.

### N3 — A check that always reports "no problem" `severity: warning · category: technical_resource`

`Execution.jl:500-518`: `check_batch_confounding` returns
`has_batch_confounding => false` unconditionally, with
`note => "Stub: real implementation would check correlation between batch and
group via chi-square or ANOVA"`.

The note is in the manifest, so the information is not lost — but the key a
reader, a UI, or an export gate will look at is `has_batch_confounding`, and it
says `false`. That is the difference between "we checked" and "we did not".

**Change (low risk, wording + one field):** add
`"checked" => false` (or `status => "not_implemented"`) and classify the notice
as `technical_resource`. Do **not** delete the check — deleting it would remove
the record that it is owed. **Belongs in:** MetaManifold.

### N4 — Two different operations are both called `rarefy` `severity: danger · category: method_limitation`

| where | what it does |
|---|---|
| `diversity.jl:55-85` (`rarefy`) | genuine without-replacement subsampling: a pool of one entry per read, partial Fisher–Yates, `depth` draws |
| `Execution.jl:1225-1234` (`normalization.method = "rarefy"`) | `prepared[:, j] = counts_after_zero[:, j] .* (min_lib / lib_sizes[j])` — a comment at line 1230 says *"For stub, rarefy by subsampling proportionally to min_lib (not exact, just scaling)"* |

Both are reachable from the product: the catalogue
(`docs/statistics/method-catalogue-v1.md:32`) lists `rarefy` under
`diversity.jl`, and `AnalysisConfig.VALID_NORMALIZATION_FOR_METHOD[NB_GLM]`
(`AnalysisConfig.jl:121`) lists `rarefy` for the analysis path too. A user who
picks `normalization.method = "rarefy"` for `nb_glm` gets **scaled non-integer
counts** handed to `MASS::glm.nb` — a count model on fractional pseudo-counts,
with none of the variance behaviour that rarefaction is chosen for.

This is precisely the substitution the rest of the repository refuses
(`Estimation.REFUSED_DISPERSION`, `ilr_basis.jl` refusal list).

**Change (low risk):** do not silently change the numbers. Either (a) refuse
`normalization.method = "rarefy"` in the analysis path with a message naming
the pipeline's `rarefy` and the fact that this path would scale rather than
subsample, or (b) implement real subsampling there and record the seed. Either
way, the DANGER banner already fires for `rarefy` + `nb_glm`
(`AnalysisConfig.jl:1353`). **Belongs in:** MetaManifold (a science decision).

### N5 — Ward squares its dissimilarities under a 1e300 sentinel `severity: info · category: technical_resource`

`ilr_basis.jl:785-787` squares `diss` in place for `ward`, then divides by
member counts in the Lance–Williams update. Values above ~1e154 would square to
`Inf`. Log-ratio variances never get there, so this is a documented guard rather
than a defect — but it is undocumented in the code. One comment's worth of
change. **Belongs in:** MetaManifold.

## Backend and environment

### B1 — CI has produced no verdict on the audited code `severity: fatal · category: dependency_environment`

```
$ gh run list --limit 10
completed  startup_failure  triage(backlog): …            CI    main                 0s
completed  startup_failure  CI                            CI    arena/01a0dd12…      0s
completed  startup_failure  feat(analysis): PhILR, SBP and balance-dendrogram ILR bases …  CI  main  0s
completed  success          PR #77                        CodeQL  refs/pull/77/head  59s
completed  failure          npm_and_yarn in /test/doi …   Dependabot Updates  main   38s

$ gh api …/actions/runs/36235457659/jobs
{"total_count":0,"jobs":[]}
```

Every `CI` run since the PR #75 merge is `startup_failure` with **zero jobs**;
GitHub's own verdict is *"This run likely failed because of a workflow file
issue."* CodeQL (a `dynamic` workflow, not repo-file-based) succeeds, which
localises the fault to the workflow dispatch rather than the runner.

This is not a code finding, but it is the first of the owner's four
preconditions and it is not met. It also means the ILR bases landed with
**no test run and no benchmark run recorded anywhere** — `bench/results/` is
git-ignored (`.gitignore:281`) and empty.

**Change:** outside this repository's diff — owner action on the workflow /
repo settings. Recorded here because it gates everything in Part 5.

### B2 — RCall holds a second copy of the prepared table, and the cleanup is on the success path only `severity: warning · category: dependency_environment`

`estimation.jl:458` `RCall.globalEnv[:est_counts] = prepared` copies the whole
`features × samples` matrix into R's global environment; Julia's `prepared` is
still live. Peak is 2× the largest object for the duration of the fit.

Worse, the release is inside the `try`:

```
478:  RCall.reval("rm(list = ls(pattern = \"^est_\")); gc()")
…
488:  catch err            # anything thrown above skips line 478
498:  finally
499:      rm(dir; recursive = true, force = true)
```

Any throw between 458 and 478 — an R error inside `glm.nb`, a `stop()` from
`R_ESTIMATION_SETUP`, an interrupt — leaves `est_counts` in R's global
environment for the life of the process, because `finally` only removes the temp
directory. The next run's `rm(list = ls(pattern = "^est_"))` will eventually
clear it, but a pipeline run holding the lock can leave it resident for hours.

**Change (low risk, three lines):** move the `rm(... est_* ...); gc()` into the
`finally`, guarded by the runtime lock. The 2× copy is inherent to RCall's
hand-off and is not fixed here — see Part 4.
**Belongs in:** MetaManifold (the cleanup) / a standalone runtime package (the
hand-off protocol).

### B3 — A missing or busy R gives four different answers, and one of them is a plausible-looking result `severity: danger · category: result_safety`

Question 9, answered: **no, not uniformly.** Four paths, four behaviours:

| path | when R/vegan is unavailable | when R is busy |
|---|---|---|
| `_alpha_significance` (`analysis.jl:605-622`) | `_not_computed(:r_unavailable, reason)` — a typed "not run" with a reason, rendered by `_significance_caption` as *"Kruskal–Wallis not run<br>…"* | caught `RBusyError` → `_not_computed(:r_busy, …)` |
| `run_nmds` (`analysis.jl:1020`) | `(fill(NaN, size(mat,1), 2), NaN)` — **a full-size coordinate matrix of NaN**, and the docstring at `1015-1018` says *"Returns NaN-filled results on failure"* | `RBusyError` propagates out of `with_r_lock` — **unhandled** |
| `run_permanova` (`analysis.jl:1057-1091`) | `nothing` — the reason (R unavailable / no usable covariate / `adonis2` threw) is **discarded** | same unhandled propagation |
| the HTTP routes (`routes/analysis.jl:771`, `802`, `860`) | 503 `r_unavailable` with a clear message; NaN → 500 `nmds_failed`; `nothing` → 500 `permanova_failed` | falls through to the framework's 500 |

The routes save the user today: `any(isnan, coords)` at `routes/analysis.jl:802`
converts the NaN matrix into a 500, and `_ensure_r()` is checked before both.
But the *library* contract is wrong — `run_nmds` hands back a matrix of the
right shape full of NaN, which is exactly the "plausible replacement result"
this repository refuses elsewhere, and `run_permanova`'s `nothing` cannot
distinguish "no covariates with ≥2 levels" (a fact about the data) from "vegan
is not installed" (a fact about the machine) from "`adonis2` failed" (a fact
about the fit).

**Change (low risk):** give both functions the `AlphaSignificance` treatment — a
return type with a status and a reason, `r_unavailable` / `r_busy` /
`no_covariates` / `fit_failed` — and catch `RBusyError` in both, as
`_alpha_significance` already does. `_ensure_r` should return the
`Provenance.RProbeStatus` (see B4) rather than a `Bool`, so the reason is
carried rather than re-derived. **Belongs in:** MetaManifold.

### B4 — The provenance inventory does not include the package the estimator requires `severity: warning · category: dependency_environment`

- `provenance.jl:303` — `const R_PACKAGES = ["dada2", "Biostrings", "ShortRead", "vegan"]`.
- `estimation.jl:675` — `suppressPackageStartupMessages(library(MASS))`. **`MASS`
  is not in `R_PACKAGES`.** It *is* pinned in `renv.lock` (81 packages, `MASS`
  present), so this is an inventory gap, not a missing pin.
- Consequence: `probe_r` can return `:ok` — *every required package present* —
  on a machine where the estimator will fail on its first line. The failure
  does land in the right place (`estimation.jl:492` → `_not_run` with
  `"the fit could not be run: …"`) but it lands *after* provenance has declared
  the environment ready.
- Separately, `Execution.jl:133` — `RAdapter(; r_packages = String["DESeq2", "edgeR"])`.
  **Neither is in `renv.lock`.** Those names are copied into the manifest
  (`Execution.jl:1419`) as `adapter_config.r_packages`, so the manifest
  asserts a dependency the pinned environment does not provide and nothing
  verifies.

**Change (low risk):** add `MASS` to `R_PACKAGES`; change `RAdapter`'s default
`r_packages` to the packages the pinned environment actually carries and that
the code actually loads (`MASS`); and have the manifest record
`r_packages_verified => true/false` with the probe status rather than a bare
list. **Belongs in:** MetaManifold.

### B5 — R-side warnings raised outside the per-feature handler never reach the manifest `severity: notice · category: technical_resource`

`estimation.jl:741-748` wraps each fit in
`withCallingHandlers(..., warning = function(w) { warns <<- …; invokeRestart("muffleWarning") })`
and folds the captured text into `est_note[i]`. That is a **capture**, not a
suppression, and it is done well: a `glm.nb` "iteration limit reached" becomes a
row-level `status = "failed"` with the message attached.

But anything R warns about *outside* that handler — `library(MASS)`,
`as.data.frame(est_cols)`, `write.csv` — goes to R's default handler, into the
process log, and nowhere near the result. Likewise `analysis.jl:1031`
(`suppressWarnings(friedman.test(...))`) and the `suppressPackageStartupMessages`
calls: targeted, justified, and invisible to the manifest.

**Change (low risk):** set `options(warn = 1)` (or collect with
`withCallingHandlers` around the whole `reval`) and append anything captured
outside a feature loop to `diagnostics["r_warnings"]`. **Belongs in:**
MetaManifold.

### B6 — `JuliaAdapter` is a shape without a body `severity: notice · category: technical_resource`

`Execution.jl:159-178` defines `JuliaAdapter(method, use_multithreading, seed,
optimizer)`. `run_analysis` uses the adapter for exactly two things: the method
match check (`Execution.jl:1580`) and `seed` (`1636`). Estimation always goes
through `Estimation.estimate_models`, i.e. R. There is no pure-Julia estimator
to serve as a reference path (question 10 — qualified in Part 2).

**Change (low risk, honesty first):** document in `run_analysis`'s docstring
that `JuliaAdapter` selects no estimator today, or make constructing one throw
until it does. Do not leave a type that looks like a backend choice and is not
one. **Belongs in:** MetaManifold.

## Warning boundary and provenance

### W1 — A run that produced no statistics is reported as safe `severity: danger · category: result_safety`

`Execution.jl:1641-1658`:

```
outcome = Estimation.estimate_models(…)          # may return :not_run
…
if outcome.status == :not_run
    push!(warnings, "Estimation not run: $(outcome.reason)")   # recorded ✓
end
exec_diagnostics = ExecutionDiagnostics(
    …
    is_dangerous = diagnostics.is_dangerous,     # ← unchanged
    banner       = diagnostics.banner)           # ← unchanged
```

The reason *is* recorded, in two places (`warnings` and
`checks["estimation"]`), which is the good half. The other half:
`is_dangerous` and `banner` are copied from the pre-estimation diagnostics, so
they do not know that nothing was fitted. And
`log_danger_banner_execution` runs **before** estimation (`Execution.jl:1619`),
so the log line for such a run is:

```
┌ Info: Execution diagnostics safe — no DANGER banner   config_id = …
```

…followed by `@info "run_analysis finished" … features=0 estimation_status=not_run`.

A result object with an empty `results` dictionary, no banner,
`is_dangerous = false`, and a log saying "safe" is the exact shape that gets
exported.

**Change (low risk, five lines):** after `estimate_models`, set
`is_dangerous |= outcome.status != :ok` and append the outcome's reason to the
banner; move `log_danger_banner_execution` after estimation. **Belongs in:**
MetaManifold. This is the single highest-value change in the audit.

### W2 — Warnings and healings are preserved; disabled analyses are not `severity: warning · category: dependency_environment`

Question 8, answered in three parts:

| recorded thing | where it lands | verdict |
|---|---|---|
| healings | `diagnostics.healings`, `provenance["diagnostics"]["healings"]` | **preserved** (counts only — see N2) |
| warnings | `diagnostics.warnings`, `provenance["diagnostics"]["warnings"]` | **preserved** |
| ILR basis warnings + DANGER reasons | `checks["ilr"]`, appended to `warnings` after `self_diagnostics` (`Execution.jl:1290-1295`, `1372-1400`) | **preserved**, with a comment explaining why the append must happen after |
| scaling warnings | `checks["scaling"]` + `append!(warnings, scaling.warnings)` (`Execution.jl:1286-1289`, with a comment about not being clobbered) | **preserved** |
| estimation `:not_run` | `checks["estimation"]` + a warning | **preserved in the record, invisible to `is_dangerous`** (W1) |
| **a disabled analysis** (vegan absent → no NMDS/PERMANOVA) | `@warn` in the log, and a 503/500 from the route | **not in any manifest** — chart routes have no manifest at all |
| **a degraded chart** (boxplot with no significance annotations) | caption text `_significance_caption` | **preserved for alpha only** |

The gap is structural: anything without a manifest cannot record that it was
degraded. Closing it is what the Run Notices model is for — a notice list that
travels with the *response*, not only with the *result object*.

### W3 — Warnings are free text, so nothing can filter, group or block on them `severity: notice · category: technical_resource`

Every warning in the system is a `String` pushed onto a `Vector{String}`:

```
push!(warnings, "Low prevalence/abundance: $(…) taxa below min_prevalence=$(…) — will be filtered")
```

There is no id, no severity, no category, no "does this block export". The
frontend cannot act on them and does not try: `DangerBanner.tsx` re-derives its
own four reasons from the config independently
(`DangerBanner.tsx:15-26`), while the backend derives a different list
(`AnalysisConfig.jl:1345-1360`). Two lists, same question, no shared source —
the drift this repository's own comments warn about when a hook and a CI check
re-implement the same rule.

**Change:** the Run Notices model (Part 3 and `docs/notices/`). Staged so that
free-text warnings keep working while notices are introduced.

### W4 — Benchmark warnings conflate three different things `severity: notice · category: technical_resource`

Question 6, answered: **no.** `bench/ilr_bases/benchmark.jl:70-77`:

```
function measure(f)
    GC.gc()
    rss0 = Sys.maxrss()
    stats = @timed f()
    return (seconds=…, allocated_bytes=stats.bytes, gc_seconds=…,
            peak_rss_growth_bytes = max(0, Int(Sys.maxrss()) - Int(rss0)))
end
```

Three separate conflations:

1. **Cumulative allocation vs peak memory.** `stats.bytes` is *total bytes
   allocated* — a churn metric. A transform that allocates and frees 8 GB while
   never holding more than 5 MB reports 8 GB. The warning text
   (`benchmark.jl:118`) says *"allocated 7.45 GiB (> 1 GiB)"*, which a reader
   will read as a memory footprint. Both numbers are printed; only one is
   described in words, and it is the misleading one.
2. **Algorithm vs setup/I/O.** For `phylogenetic` and
   `sequential_binary_partition`, `ilr_transform` reads the file from disk,
   SHA-256s it, copies it into a `String`, and parses it — all inside the
   measured closure (the docstring says so: *"file read and SHA-256 included"*).
   The file *write* is outside (`workloads()` at line 79-99 runs before
   `measure`). So the SBP number is mostly I/O and text handling, and the
   balance computation is a rounding error inside it. Conversely the `clr` and
   `ilr default Helmert` workloads go through `prepare_analysis_table`, which
   includes config validation, two prevalence/abundance passes, diagnostics and
   manifest construction — so those two are full-pipeline numbers compared
   against engine-only numbers in the same table.
3. **Peak RSS measured with a monotone high-water mark.** `Sys.maxrss()` never
   decreases. After the first workload grows the heap, every later workload that
   reuses it reports `Δ = 0`, which reads as "used no memory" when it means "used
   no *new* memory".

**Change (low risk, benchmark-only):** decompose `measure` into labelled phases
(`source_read`, `source_hash`, `basis_build`, `balances`, `checks`) — the
functions are already separate, so this is wrapping, not refactoring; report
`allocated_bytes` and `peak_rss_growth_bytes` as distinct named quantities with
distinct thresholds; and word the warnings as *"allocated X (churn) / grew peak
RSS by Y (footprint)"*. **Belongs in:** MetaManifold (bench scripts) — though a
shared phase-timer is a good candidate for the standalone runtime package
later.

### W5 — No global suppression exists `severity: n/a · category: n/a`

Confirmed by search, since the brief asks for it. The only places logging is
intercepted:

| site | what it does |
|---|---|
| `bench/ilr_bases/benchmark.jl:108`, `regression_gate.jl:62`, `test/unit/test_provenance.jl:262`, `src/doi/Zenodo.jl:33` | `NullLogger` around a measured or network block — benchmark noise, not a warning policy |
| `src/server/server.jl:216` | `_SuppressEpipe(global_logger())` — filters broken-pipe errors from HTTP handlers only |
| `estimation.jl:741-748` | R `muffleWarning` around a single fit — **captures** the warning text into `est_note[i]`, so nothing is lost |
| `analysis.jl:1031`, `estimation.jl:674` | `suppressWarnings(friedman.test(...))`, `suppressPackageStartupMessages` — targeted to one call each |

None is a global suppression. `run_nmds`'s `tryCatch(..., error = function(e) NULL)`
at `analysis.jl:1026-1034` is the closest thing to a silent swallow — it discards
the error *message* and returns NULL, which is then turned into NaN coordinates
(see B3).

### W6 — Two warnings that repeat, and one check computed twice `severity: info · category: technical_resource`

- `_ensure_r()` (`analysis.jl:1000-1008`) logs `@warn "R/vegan not available -
  NMDS and PERMANOVA disabled"` **on every call** while vegan is absent, because
  `_r_loaded[]` stays `false`. Each call also takes the R lock (10 s timeout) and
  runs a failing `library(vegan)`. Harmless once, noisy under a polling UI.
- `check_prevalence_abundance` runs twice per analysis: inline at
  `Execution.jl:848-856` (on the filtered matrix, to build `keep_taxa`) and again
  inside `self_diagnostics` (`Execution.jl:545`, deliberately on `raw_counts`, to
  report how many taxa *will be* filtered). Both allocate a row copy per taxon
  (`counts[i, :]`); `@view` or `sum(>(0), @view counts[i, :])` removes D
  allocations per call.

**Change (low risk):** memoise `_ensure_r`'s negative result with the reason;
use views in `check_prevalence_abundance` and the two `check_all_zero_*` loops.
**Belongs in:** MetaManifold.

---

# Part 2 — The ten questions, answered

1. **ILR/PhILR/SBP/balance-dendrogram memory behaviour.** The engine is O(D·n)
   and builds no dense basis, as documented — verified: `tree_balances`
   (`ilr_basis.jl:1078-1105`) is one post-order pass per sample over three
   `2D−1` vectors, and the SBP support sets are keyed by `(size, hash)` rather
   than by index vectors precisely to avoid O(D²) keys. The peak costs are
   *around* the engine, not in it: the `diss` copy (**M2**), the SBP source
   materialisation (**M5**), and the `clr_table` (**M3**). Separately, the
   `default` basis that most users actually get is O(D²·n) in allocation
   (**M1**).

2. **Are dense intermediates materialised unnecessarily?** Yes, four of them:
   `clr_table` (**M3**), the `diss` copy (**M2**), the `String(copy(bytes))`
   (**M5**), and the DataFrame→Matrix duplication in `filtered_counts` (**M6**).
   Plus three redundant full-table copies in `prepare_analysis_table` (**M4**).
   None is needed.

3. **Would block processing, views, preallocation or streaming help?** Yes, in
   this order of value: views in the row-wise diagnostic loops (**W6**);
   preallocation for the CLR/ILR column vectors (**M3**); streaming the SBP
   straight from an `IO` into `W` with a chunked hash (**M5**); block
   (column-block) processing only if the profiling says the `variation_condensed`
   kernel matters (**M8**). Nothing here needs a new library.

4. **Does RCall duplicate large data structures?** Yes, in two ways. Inherently:
   `RCall.globalEnv[:est_counts] = prepared` (`estimation.jl:458`) copies the
   whole table into R while Julia keeps its own, so peak is 2× the largest
   object for the duration of the fit. Avoidably: the `rm(...); gc()` that
   releases it is inside the `try`, so any throw leaks it for the process
   lifetime (**B2**). The per-feature loop then extracts one row at a time
   (`y <- as.numeric(est_counts[i, ])`), which is the right shape.

5. **Can DuckDB or on-disk intermediates reduce peak memory?** DuckDB is
   *already* the on-disk IR for the analysis input — `filtered_counts` runs SQL
   against `results.duckdb`. The peak problem is the hand-off, not the storage:
   the query result is materialised as a `DataFrame` and then copied into a
   transposed `Matrix` cell by cell (**M6**). Fixing that is worth more than
   routing the transform through SQL. For the transform path specifically, the
   dominant peak is the `D(D−1)/2` variation vector, which is intrinsically
   dense and sample-major-unfriendly; spilling it to disk is possible but would
   be a different algorithm with different numerics, and the honest first move
   is to halve it (**M2**) and then to reconsider the threshold.

6. **Do benchmark warnings distinguish algorithmic memory from setup/I/O
   memory?** No (**W4**). Three conflations: cumulative allocation vs peak
   footprint; file read + hash + parse inside the measured closure vs file write
   outside; and `Sys.maxrss()` as a monotone high-water mark.

7. **Does self-healing change data or only metadata?** **Data.** NaN/Inf
   healing rewrites values in `prepared`; `IMPUTE` writes `epsilon` into
   `counts` before the transform; `DROP` removes rows/columns. What is recorded
   is a sentence with a count, never the positions — so healed values are
   indistinguishable from measured ones (**N1**, **N2**).

8. **Is every healing, fallback, disabled analysis and warning preserved in the
   result manifest?** Warnings, healings, scaling warnings, ILR warnings and
   DANGER reasons, and the estimation `:not_run` reason: **all preserved**.
   Disabled analyses (vegan missing → no NMDS/PERMANOVA) and degraded charts:
   **not preserved anywhere**, because the routes that produce them have no
   manifest (**W2**). And `is_dangerous` does not learn about a `:not_run`
   estimation (**W1**).

9. **Does a missing R package cause a clear "not run" rather than a plausible
   replacement?** For the *estimator*, yes — `Estimation` returns
   `EstimationOutcome(:not_run, reason, …)` with an empty `results`
   (`estimation.jl:263-275`, `488-495`), and there is deliberately no fourth
   status meaning "something plausible was produced". For the *ordination and
   significance* layer, no — four different behaviours, one of which is a
   NaN-filled coordinate matrix (**B3**). And the provenance probe does not
   check `MASS`, the package the estimator loads (**B4**).

10. **Can the pure Julia implementation serve as a reference path for future
    accelerated implementations?** **For the transform, yes — and it already
    does.** `ILRBasis` is pure Julia, calls no R, and is validated three ways:
    Agda proofs (`proofs/agda`, cited by theorem name in `ilr_basis.jl:16-27`),
    an independent Base-only port of `philr`/`compositions`/`robCompositions`
    (`test/fixtures/ilr/ilr_reference.jl`, Base Julia only "so that nothing here
    can drift with a dependency", agreeing to 1e-10 on committed fixtures), and
    the `comb-is-helmert` test tying it to the historical default at 1e-12.
    That is exactly the reference a future kernel needs: an O(D·n) algorithm, a
    tolerance, and an oracle.
    **For the estimator, no.** `Estimation` is R-only and `JuliaAdapter`
    (`Execution.jl:159-178`) is a declared type with no implementation behind it
    (**B6**). A reference path for a future accelerated *estimator* does not
    exist yet and should be a precondition of building one.

---

# Part 3 — Classification of every warning found

Category key: **T** technical/resource · **D** dependency/environment ·
**Q** data quality · **M** method limitation · **R** result safety ·
**F** fatal. Severity: info · notice · warning · danger · fatal.

## Harmless engineering noise

| id | site | today | category | sev | action |
|---|---|---|---|---|---|
| `ilr.part_weights_unused_kinds` | `ilr_basis.jl:976` | silent | T | info | compute only the selected kind |
| `exec.prepared_copy_unused` | `Execution.jl:1059` | silent | T | info | move the copy into `clr`/`rarefy` |
| `exec.counts_copy_unused` | `Execution.jl:840` | silent | T | info | drop |
| `bench.phase_undecomposed` | `bench/ilr_bases/benchmark.jl:70` | `::warning` | T | notice | phase the measurement |
| `analysis.r_unavailable_relogged` | `analysis.jl:1000` | `@warn` per call | D | notice | memoise |

## Genuine environment problems

| id | site | today | category | sev | action |
|---|---|---|---|---|---|
| `ci.no_verdict` | `.github/workflows/ci.yml` | none | D | **fatal** | owner action; blocks the whole gate |
| `r.package_not_probed` (`MASS`) | `provenance.jl:303` vs `estimation.jl:675` | none | D | warning | add `MASS` to `R_PACKAGES` |
| `r.adapter_packages_unpinned` | `Execution.jl:133` | manifest claims them | D | warning | default to pinned packages; mark verified |
| `r.est_counts_leaked_on_error` | `estimation.jl:478` | none | D | warning | clean up in `finally` |
| `r.warning_outside_handler` | `estimation.jl:741` | process log only | T | notice | collect into diagnostics |
| `r.unavailable` | `analysis.jl:1005` | `@warn` | D | warning | return `RProbeStatus`, record as a notice |
| `r.busy` | `analysis.jl:617`, `1023` | `@warn` / unhandled | D | warning | catch in both ordination paths |

## Scientific interpretation warnings

| id | site | today | category | sev | notes |
|---|---|---|---|---|---|
| `ilr.dendrogram_data_derived` | `ilr_basis.jl:1265` | in `checks` only | M | warning | *"the basis is chosen from the same table that is then tested"* — already worded well; it needs a notice so the UI shows it |
| `ilr.sbp_p_hacking` | `ilr_basis.jl:1245` | DANGER banner | R | **danger** | correct as-is |
| `ilr.balance_weights_not_isometric` | `ilr_basis.jl:1291` | in `checks` only | M | warning | *"effect sizes change, per-balance test statistics do not"* |
| `ilr.pruned_tips` | `ilr_basis.jl:1225` | warning string | M | warning | |
| `ilr.zero_length_tip_edges_replaced` | `ilr_basis.jl:1271` | warning string | M | warning | |
| `exec.rarefy_discouraged` | `Execution.jl:1220` | `@warn` + DANGER | M | danger | see N4 — the method is not rarefaction |
| `config.pseudocount_unusual` (≥1, <0.1) | `AnalysisConfig.jl:182`, `401-404` | `@warn` | M | warning | |
| `config.epsilon_extreme` | `AnalysisConfig.jl:212`, `412-415` | `@warn` | M | warning | |
| `config.css_quantile_below_median` | `AnalysisConfig.jl:229` | `@warn` | M | warning | |
| `config.tmm_trims_zero` | `AnalysisConfig.jl:240` | `@warn` | M | warning | |
| `config.parameter_inert` (`css_quantile` / `tmm_*` with a method that ignores them) | `AnalysisConfig.jl:249-252` | `@warn` | M | warning | a recorded parameter that does nothing |
| `scaling.*` (CSS quantile, TMM undefined, RLE zero features) | `scaling.jl:214`, `342`, `348`, `410`, `416` | warning strings | Q | warning | already well worded |
| `estimation.boundary_theta` / `separation` | `estimation.jl:786-799` | row `status = "boundary"` | M | warning | correctly not counted as `ok` |

## Possible correctness risks

| id | site | risk | category | sev |
|---|---|---|---|---|
| `estimation.not_run_reported_safe` | `Execution.jl:1641-1658` | zero results, `is_dangerous = false`, log says "safe" | R | **danger** |
| `nmds.nan_filled_result` | `analysis.jl:1020` | a full-size NaN coordinate matrix is a plausible-looking result | R | **danger** |
| `permanova.reason_discarded` | `analysis.jl:1086-1091` | `nothing` conflates three distinct causes | R | warning |
| `exec.batch_confounding_stub` | `Execution.jl:500` | renders as "no confounding found" | T | warning |
| `heal.nan_inf_scale_mixed` | `Execution.jl:600-616` | writes out-of-scale values into the fitted table | R | **danger** |
| `heal.positions_unrecorded` | `Execution.jl:706` | healed cells indistinguishable from measured ones | R | warning |
| `exec.rarefy_is_scaling` | `Execution.jl:1225-1234` | fractional counts into a count model, under the name "rarefy" | M | **danger** |
| `ilr.dendrogram_threshold_half` | `ilr_basis.jl:1250-1252` | refuses at a reported 2 GiB that is really 4 GiB | T | warning |
| `ilr.helmert_quadratic` | `Execution.jl:1188` | 8 GB of churn at the benchmark's largest size; not a wrong answer, but a wall | T | notice |

---

# Part 4 — Where each change belongs

The brief asks for three buckets. The test I applied: **does the change encode
a scientific decision about this product's results, or is it infrastructure
with its own test matrix and no semantics?**

## MetaManifold (this repository)

Everything below is a change to what the product *says* or *records*, needs no
new dependency, and is low-risk by construction.

1. **Classification and provenance** — add `id`, `category`, `severity`,
   `action_required`, `result_usable`, `blocks_export`, `log_ref` to every
   notice (Part 3's tables are the inventory to start from).
2. **W1** — `is_dangerous` and the banner must learn that estimation did not
   run; `log_danger_banner_execution` moves after estimation.
3. **B3** — give `run_nmds` and `run_permanova` typed not-run outcomes, and
   catch `RBusyError` in both, as `_alpha_significance` already does.
4. **N3** — mark `check_batch_confounding` as not-implemented rather than
   negative.
5. **B4** — add `MASS` to `R_PACKAGES`; make `RAdapter`'s default package list
   the pinned one; record verification status in the manifest.
6. **B2** — release R's `est_*` objects in a `finally`.
7. **Memory, in this order:** M2 (the `diss` copy and the threshold constant),
   M3 (`clr_table`), M4 (three copies), M6 (DuckDB hand-off), M8 (locality),
   M7 (unused weight kinds). Each is a few lines and none changes a number.
8. **W4** — decompose `measure` into phases; report churn and footprint as
   separate named quantities.
9. **Wording** — the benchmark warning ("allocated" → "allocated churn" vs
   "grew peak RSS"), the dendrogram refusal message (state the true peak), and
   `run_nmds`'s docstring (it currently promises NaN-filled results).
10. **Documentation** — this audit; the Run Notices model
    (`docs/notices/`); a note in `docs/statistics/method-catalogue-v1.md`
    distinguishing the two `rarefy` implementations.

## A standalone runtime package

These are infrastructure. They have large, boring test matrices (allocation
bounds, streaming edge cases, concurrency) that would drown a statistics suite,
and no scientific semantics of their own.

1. **The R hand-off protocol.** The 2× copy at `estimation.jl:458` is inherent
   to RCall. A zero-copy or shared-memory hand-off (Arrow C data interface, or
   a shared array for the duration of the fit) belongs here, together with the
   session-pool problem implied by `r_runtime.jl` (one embedded interpreter for
   the whole process, serialised by a single `ReentrantLock`).
2. **Streaming readers with a bounded-allocation contract** — Newick and SBP
   CSV consumed from an `IO` with a documented peak, plus chunked SHA-256. This
   is M5(b) and M5(c), and it wants a property-based suite
   ("peak ≤ c·D for any input"), not a statistics test.
3. **The notice runtime** — the `RunNotice` type, the collector, the
   JSON/Nickel schema, log correlation ids, and the "does anything here block
   export" fold. MetaManifold would depend on it and supply the catalogue.
4. **A spill-to-disk variation matrix / blockwise `hclust`**, if the profiling
   says the dendrogram path needs it. Different algorithm, same numerics to
   prove — that is a package-sized job.
5. **The phase-timing harness** the decomposed benchmark wants.

## Axiom or another upstream library

Nothing, today. The four preconditions are not met (B1 alone stops it), and the
brief rules it out for this task. For the record, the kernels that would
eventually belong there, and *why there*:

1. **The balance kernel** — `b_n = c_n · (mean_p over N₊ − mean_p over N₋)`
   over a post-order traversal, parameterised by element type and device. It is
   already proved (Agda), already O(D·n), and already has an oracle; the only
   thing a faster backend changes is the constant. That is the ideal
   first kernel — and the reason to do M1 first, so the default path and the
   reference path become the same code.
2. **The log/mean reductions with configurable accumulation order** — the
   parity tolerance for any accelerated kernel is meaningless unless the
   summation order is specified. This belongs with the kernel, not in
   MetaManifold.
3. **Arbitrary-precision / exact kernels** for `NumericPolicy`'s
   `:exact_counts` and `:high_precision` modes, where the current code already
   names the seam (`with_precision`, `ResourceLimitError`).
4. **A batched `hclust`** with the same tie-breaking as R's `hclust.f` — but
   only after the tie behaviour is pinned down by tests, which is a MetaManifold
   job first.

The general rule: **MetaManifold owns the conditions, the refusal policy, the
tolerances and the oracles; the kernel library owns the constant factor.** The
repository already works this way — `ilr_basis.jl:7-9` says the conditions
document outranks the implementation, which is exactly the boundary.

---

# Part 5 — Recommended order of work

Everything before the gate is additive, changes no number, and can land on a
branch while CI is broken.

### Stage 0 — unblock the evidence (owner action, not code)

1. Fix the `startup_failure` (B1) and get one green `CI` run on `main`.
2. Run `julia --project=. bench/ilr_bases/benchmark.jl` and commit the JSON —
   or start publishing it as an artefact — so the DERIVED numbers in this audit
   become measurements. Until then, M1/M2/M5 are well-founded predictions, not
   findings of fact.

### Stage 1 — classification and provenance (no behaviour change)

3. Add the `RunNotice` type and the catalogue (Part 3's three tables, migrated
   into `docs/notices/catalogue.md`); emit notices alongside the existing free
   text; leave `Vector{String}` intact so nothing downstream breaks.
4. W1 — `is_dangerous`, the banner, and the position of
   `log_danger_banner_execution`.
5. N3 — the not-implemented check. B4 — `MASS` and the adapter defaults.
6. B2 — R cleanup in `finally`. B3 — typed not-run outcomes for NMDS/PERMANOVA.
7. N2 — record what healing changed.

### Stage 2 — memory, cheapest first

8. M2 (two lines: the `diss` copy, and the threshold constant).
9. M3, M4, M6, then M8 and M7.
10. Re-run the benchmark with the new phase decomposition (W4) and publish the
    before/after.

### Stage 3 — decisions that need the owner

11. N1 (how to heal, or refuse) and N4 (which `rarefy` is which). Both change
    results that saved analyses were computed from, so both want the treatment
    `docs/statistics/behaviour-change-zero-depth-samples.md` gives its own
    change: a published note, before the code.

### Then, and only then, the gate

The brief's four preconditions, with what "met" looks like:

| precondition | how it is met |
|---|---|
| current tests pass | green `CI` on `main` (B1) — **not met** |
| current CI is operational | a workflow that produces jobs — **not met** |
| reference outputs are captured | `bench/results/ilr_bases_results.json` committed or published as an artefact, plus the fixture agreement already enforced by `ilr_reference.jl` — **partially met** (fixtures yes, benchmarks no) |
| parity criteria are defined | *max abs difference < 1e-12 between the default Helmert path and `tree_balances(comb_tree(taxa), X, 1, 1)`, on the existing `comb-is-helmert` workloads and the three benchmark sizes, with the prepared-table hash unchanged for the `clr` path* — **not yet written down** |

When all four hold, M1 is the first change to make: it makes the default path
and the reference path the same code, which is the precondition for ever
asking a kernel library to make it faster.

---

## Appendix A — Derived arithmetic

`Float64 = 8 bytes`. All figures DERIVED from source, not measured.

**A1. Default Helmert allocation** (`Execution.jl:1188`).
Per sample: Σ_{i=1}^{D−1} 8i = 4·D·(D−1) bytes. Per run: × n.

| D × n | bytes |
|---|---|
| 1 000 × 20 | 8.0e7 (80 MB) |
| 1 500 × 40 | 3.6e8 (360 MB) |
| 10 000 × 20 | 8.0e9 (8.0 GB) |

**A2. Balance dendrogram peak** (`ilr_basis.jl:773`, `803`).
`tau = 8·D(D−1)/2`; `diss` = the same again; `C = 8·n·D`.
Reported `need` = `tau` only, so peak/reported = 2.00 (plus `C`).

| D | reported | peak | first D that warns | first D refused |
|---|---|---|---|---|
| 10 000 | 381 MiB | 764 MiB | — | — |
| 16 385 | 1.00 GiB | 2.00 GiB | ✓ | — |
| 23 171 | 2.00 GiB | 4.00 GiB | — | ✓ |

**A3. SBP source** (`ilr_basis.jl:1228-1231`), ~2 bytes per CSV cell.
file text = 2·D(D−1); `String(copy(bytes))` = the same; `W` (Int8) = D(D−1).

| D | text | copy | W | peak | cells |
|---|---|---|---|---|---|
| 10 000 | 191 MiB | 191 MiB | 95 MiB | 477 MiB | 1.0e8 |

**A4. Count-table copies** (`Execution.jl:840`, `858`, `868`, `1059`, `1154`),
input 20 000 taxa → 10 000 retained, 200 samples. One retained copy = 15.3 MiB;
input = 30.5 MiB. Current peak ≈ 107 MiB; after removing the dead copy, the
composable slice and `clr_table` ≈ 76 MiB.

**A5. `filtered_counts`** (`analysis.jl:112-123`), 10 000 features × 500
samples: a 40 MiB `DataFrame` plus a 40 MiB `Matrix`, filled by 5 000 000
`DataFrameRow` name lookups.

## Appendix B — What this audit did not verify

| claim | how to verify |
|---|---|
| M1's 8 GB, M2's 2×, M5's 477 MiB | `julia --project=. bench/ilr_bases/benchmark.jl ILR_BENCH_TAXA=10000` (already reports `allocated_bytes` and ΔpeakRSS) |
| that `hclust_r`'s `diss` copy is not elided | `@allocated` on `hclust_r(tau, D, "average")` vs `2·sizeof(tau)` |
| that the R copy of `est_counts` doubles peak | `Sys.maxrss()` immediately before and after `estimation.jl:458` |
| whether `MASS` is present in the CI image | `Rscript -e 'packageVersion("MASS")'`, or extend `probe_r` (B4) and read the manifest |
| the `startup_failure` root cause | `gh api /repos/hyperpolymath/MetaManifold-WebUI/actions/runs/<id>` — zero jobs, and GitHub's own "workflow file issue" verdict |
