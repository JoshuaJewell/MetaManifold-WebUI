<!--
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Run Notices — catalogue

**Status: proposed, 2026-09-26.** Every entry is `status: proposed` until the
type lands; none is emitted today. The field model is in
`docs/notices/README.md`; the audit that produced this inventory is
`docs/audit/2026-09-26-memory-numerics-warning-audit.md` (referenced below as
§M1, §N1, §B3, §W1 …).

**Rules for this file** (enforced by `test/unit/test_notices.jl`):

- an `id` is immutable; to change what a notice means, add a new `id` and mark
  the old one `status: retired`;
- every `id` here is emitted somewhere in `src/`, and every `id` emitted in
  `src/` appears here;
- values never appear in an `id` — they go in `data`;
- `usable = no` ⟹ `blocks = yes`; `sev = fatal` ⟹ `usable = no`.

**Category key.** `T` technical/resource · `D` dependency/environment ·
`Q` data quality · `M` method limitation · `R` result safety · `F` fatal.
**Severity.** `info` · `notice` · `warning` · `danger` · `fatal`.
`banner` = rendered on every figure (derived: `danger`/`fatal`, or `M` at
`warning`). `manifest` = written into the run manifest.

---

## `config.*` — declared configuration

Raised by `AnalysisConfig` when the *declared* analysis is legal but says
something the reader will not expect.

| id | cat | sev | usable / blocks | banner | summary | action | site · audit |
|---|---|---|---|---|---|---|---|
| `config.pseudocount_unusual_high` | M | warning | yes / no | yes | `pseudocount = {pseudocount} is ≥ 1; 0.5–0.65 is typical for CLR/ILR.` | Lower `normalization.pseudocount`, or record the choice. | `AnalysisConfig.jl:182, 401` |
| `config.pseudocount_very_small` | M | warning | yes / no | yes | `pseudocount = {pseudocount} is < 0.1; zeros become extreme log-ratios.` | Raise it, or use `multiplicative_replacement`. | `AnalysisConfig.jl:404` |
| `config.epsilon_large` | M | warning | yes / no | yes | `epsilon = {epsilon} is > 1e-3; it will move zero handling and every log transform.` | Lower it. | `AnalysisConfig.jl:212, 412` |
| `config.epsilon_tiny` | M | warning | yes / no | yes | `epsilon = {epsilon} is < 1e-12; it can underflow.` | Raise it. | `AnalysisConfig.jl:415` |
| `config.css_quantile_below_median` | M | warning | yes / no | yes | `css_quantile = {q} is below the median, so the cumulative sum covers less than half of each sample.` | Use 0.5–0.75 (Paulson et al. 2013). | `AnalysisConfig.jl:229` |
| `config.tmm_trims_zero` | M | warning | yes / no | yes | `tmm_log_ratio_trim = 0 and tmm_sum_trim = 0: the trimmed mean is an untrimmed mean.` | Restore trimming; that robustness is why TMM is used. | `AnalysisConfig.jl:240` |
| `config.parameter_inert` | M | warning | yes / no | yes | `{parameter} is recorded but has no effect for normalization.method = '{method}'.` | Remove it, or change the method. | `AnalysisConfig.jl:249, 252` |
| `config.max_features_small` | M | warning | yes / no | yes | `max_features = {n}; only {n} features will be tested.` | Raise it, or say why the bound is deliberate. | `AnalysisConfig.jl:437` |
| `config.min_samples_per_group_low` | M | danger | yes / no | yes | `min_samples_per_group = {n} < 3; variance estimation is unstable.` | Increase the group size. | `AnalysisConfig.jl:446` |
| `config.outcome_column_undeclared` | M | warning | yes / no | no | `outcome_column '{col}' is not listed in metadata_columns.` | List it. | `AnalysisConfig.jl:647` |
| `config.outcome_column_inert` | M | warning | yes / no | no | `outcome_column is set but method '{method}' is not logistic; it will be ignored.` | Remove it, or use `logistic`. | `AnalysisConfig.jl:651` |
| `config.rarefy_with_nb_glm` | M | danger | yes / no | yes | `rarefy + nb_glm: rarefaction discards data and NB_GLM already models depth through size_factors.` | Use `size_factors` or `tss`. | `AnalysisConfig.jl:1353` |
| `config.zero_handling_refuse` | R | danger | yes / yes | yes | `zero_policy = 'refuse' with zeros present.` | Use `pseudocount` or `multiplicative_replacement`. | `AnalysisConfig.jl:184` |
| `config.correction_disabled` | R | danger | yes / yes | yes | `BH correction is disabled; false discoveries are not controlled.` | Re-enable, or acknowledge with `DANGER_ACK_TOKEN`. | `AnalysisConfig.jl:1353` |
| `config.dangerous` | R | danger | yes / no | yes | `{n} configuration overrides are in force; the DANGER banner applies.` | Review each reason before publishing. | `AnalysisConfig.jl:892` |

## `exec.*` — preparation, diagnostics, healing

| id | cat | sev | usable / blocks | banner | summary | action | site · audit |
|---|---|---|---|---|---|---|---|
| `exec.heal.nan_inf` | R | danger | **no / yes** | yes | `{nan_count} NaN and {inf_count} Inf in the prepared table were replaced with {rule}.` | Choose `drop_policy = refuse`, or re-run with a transform that does not produce them. | `Execution.jl:600, 1305` · §N1, §N2 |
| `exec.heal.nan_inf.value_scale_mixed` | R | danger | **no / yes** | yes | `Healing replaced NaN with {epsilon} but ±Inf with log(1/{epsilon}); those are different scales for one table.` | Pick one scale and record it. | `Execution.jl:600-616` · §N1 |
| `exec.heal.all_zero_samples_dropped` | Q | warning | yes / no | yes | `{n} samples with no reads at all were dropped before the transform.` | Re-run; their relative abundances do not exist. | `Execution.jl:654` · §N2 |
| `exec.heal.all_zero_samples_imputed` | M | warning | yes / no | yes | `{n} samples with no reads at all were imputed with {epsilon} before the transform.` | Prefer `drop`; imputation invents a uniform distribution. | `Execution.jl:654` · §N2 |
| `exec.heal.all_zero_taxa_dropped` | Q | warning | yes / no | yes | `{n} all-zero taxa were dropped after the transform.` | Note that balances and CLR centres were computed with them present. | `Execution.jl:675, 1320` · §N2 |
| `exec.heal.all_zero_taxa_imputed` | M | warning | yes / no | yes | `{n} all-zero taxa were imputed with {epsilon}.` | Prefer `drop`. | `Execution.jl:675` · §N2 |
| `exec.heal.positions_unrecorded` | R | warning | yes / no | yes | `{n} values were rewritten but only the count was recorded, so healed cells cannot be told from measured ones.` | Re-run once positions are recorded. | `Execution.jl:706` · §N2 |
| `exec.diag.zero_variance` | Q | warning | yes / no | no | `{n_taxa} taxa and {n_samples} samples have zero variance; singularities are likely in LM/GLM.` | Inspect them. | `Execution.jl:556` |
| `exec.diag.library_size_outliers` | Q | warning | yes / no | no | `{n} samples are more than 3 sd from the mean library size.` | Check for contamination or a failed run. | `Execution.jl:563` |
| `exec.diag.low_prevalence_abundance` | Q | warning | yes / no | no | `{n_prev} taxa below min_prevalence and {n_abund} below min_abundance were filtered.` | Relax the thresholds if the biology is rare. | `Execution.jl:568` |
| `exec.diag.batch_confounding_not_implemented` | T | warning | yes / no | no | **Batch confounding was not checked; this check is not implemented.** | Treat confounding as unknown. | `Execution.jl:500` · §N3 |
| `exec.filter.max_features` | M | warning | yes / no | yes | `Filtered to the {n} most abundant taxa; the rest were not tested.` | Raise `max_features`. | `Execution.jl:870` |
| `exec.epistemic.no_avec_fibre` | Q | warning | yes / no | yes | `No taxon has avec_fibre = true, so epistemic filtering would remove everything; nothing was filtered.` | Check the taxon metadata. | `Execution.jl:883` |
| `exec.rarefy_is_scaling` | M | danger | yes / no | yes | `normalization.method = 'rarefy' scaled each sample by min_lib/lib; it did not subsample reads.` | Use the pipeline's `rarefy` (real subsampling) or a depth offset. | `Execution.jl:1225-1234` · §N4 |
| `exec.estimation_not_run` | R | danger | **no / yes** | yes | `No statistics were produced: {reason}.` | Fix the cause named in the reason and re-run. | `Execution.jl:1641-1658` · §W1 |
| `exec.memory.allocation_churn` | T | notice | yes / no | no | `Preparation allocated {bytes} across {phases} phases.` | None; recorded for provenance. | new · §M1–M4 |
| `exec.memory.peak_rss` | T | notice | yes / no | no | `Peak process memory grew by {bytes} during preparation.` | None; recorded for provenance. | new · §M2, §M5 |

## `ilr.*` — ILR bases

| id | cat | sev | usable / blocks | banner | summary | action | site · audit |
|---|---|---|---|---|---|---|---|
| `ilr.dendrogram.variation_matrix_size` | T | warning | yes / no | no | `The variation matrix for {taxa} taxa takes {bytes}; the peak including hclust's working copy is {peak_bytes}.` | Filter taxa, or use a phylogenetic or SBP basis. | `ilr_basis.jl:1250-1252` · §M2 |
| `ilr.dendrogram.data_derived` | M | warning | yes / no | yes | `The basis was chosen from the same table that is then tested ({method} on the variation matrix).` | Pre-register the method, or use a basis chosen without reference to this table. | `ilr_basis.jl:1265` |
| `ilr.sbp.p_hacking_guard` | R | danger | yes / yes | yes | `{n} distinct SBP matrices have been tried in this project (more than {threshold}).` | Pre-register the SBP and report every partition tried. | `ilr_basis.jl:1245` |
| `ilr.phylo.pruned_tips` | M | warning | yes / no | no | `{n} tree tips that are not retained taxa were pruned (ape::keep.tip semantics) before the basis was built.` | Confirm the pruning is the one you meant. | `ilr_basis.jl:1225` |
| `ilr.balance_weights.zero_length_tip_edges` | M | warning | yes / no | no | `{n} zero-length tip edge(s) were replaced by the smallest non-zero edge length ({len}), as philr does.` | None, if that is the intended convention. | `ilr_basis.jl:1271` |
| `ilr.balance_weights.not_isometric` | M | warning | yes / no | yes | `ilr_balance_weights = '{kind}': effect sizes change, per-balance test statistics do not, and the coordinates are not isometric.` | Do not compare balance magnitudes across balances. | `ilr_basis.jl:1291` |
| `ilr.source.materialised` | T | notice | yes / no | no | `The {kind} source was held as {bytes} of file bytes plus a String copy plus a {w_bytes} matrix.` | None; recorded for provenance. | `ilr_basis.jl:1218-1231` · §M5 |
| `ilr.basis.refused` | F | fatal | **no / yes** | yes | `The {basis} basis was refused: {reason}.` | Fix the input named in the reason. | `ilr_basis.jl:71` |

## `scaling.*` — depth offsets

| id | cat | sev | usable / blocks | banner | summary | action | site · audit |
|---|---|---|---|---|---|---|---|
| `scaling.css.quantile_below_median` | M | warning | yes / no | yes | `css_quantile = {q}: the cumulative sum covers less than half of each sample.` | Use 0.5–0.75. | `scaling.jl:214` |
| `scaling.tmm.undefined` | Q | warning | yes / no | no | `The weighted trimmed mean was undefined for {samples}.` | Inspect those samples. | `scaling.jl:342` |
| `scaling.tmm.few_log_ratios` | Q | warning | yes / no | no | `Fewer than two log-ratios survived trimming for at least one sample.` | Relax the trims or check the sample. | `scaling.jl:348` |
| `scaling.rle.zero_features` | Q | warning | yes / no | no | `{n} of {total} features contain a zero and cannot enter a median-of-ratios factor.` | Note that they are excluded from the size factors. | `scaling.jl:410` |
| `scaling.rle.few_positive` | Q | warning | yes / no | no | `Only {n} features are positive in every sample.` | Note the narrow basis for the size factors. | `scaling.jl:416` |

## `estimation.*` — the fit

| id | cat | sev | usable / blocks | banner | summary | action | site · audit |
|---|---|---|---|---|---|---|---|
| `estimation.not_run` | R | danger | **no / yes** | yes | `No model was fitted: {reason}.` | Fix the cause named in the reason. | `estimation.jl:263, 488` |
| `estimation.feature_failed` | M | notice | yes / no | no | `{n} of {total} features produced no fit; each row says why and is excluded from the BH family.` | None, unless the rate is high. | `estimation.jl:~530` |
| `estimation.feature_boundary` | M | warning | yes / no | yes | `{n} features are at a boundary ({reasons}); their standard errors are not trustworthy.` | Do not interpret those rows. | `estimation.jl:786-799` |
| `estimation.dispersion_refused` | F | fatal | **no / yes** | yes | `dispersion_method '{method}' has no implementation here: {reason}.` | Use `parametric`. | `estimation.jl:70-75` |
| `estimation.exclusion_rate_high` | R | danger | yes / yes | yes | `{pct}% of features produced no fit; that is a finding about the data, not a detail.` | Investigate before publishing. | new |

## `r.*` — the embedded R runtime

| id | cat | sev | usable / blocks | banner | summary | action | site · audit |
|---|---|---|---|---|---|---|---|
| `r.unavailable` | D | error | **no / yes** | yes | `The embedded R runtime is not available, so nothing that needs it ran.` | Install R. | `analysis.jl:1005` · §B3 |
| `r.package_missing` | D | error | **no / yes** | yes | `R is reachable but {packages} are not installed.` | Restore the renv library. | `provenance.jl:358` · §B4 |
| `r.package_not_probed` | D | warning | yes / no | no | `{package} is required by the estimator but is not in the probed package list, so provenance cannot see it missing.` | Add it to `R_PACKAGES`. | `provenance.jl:303` · §B4 |
| `r.runtime_busy` | D | warning | **no / yes** | yes | `The R runtime was busy for {waited}s (a pipeline run holds it); nothing was computed.` | Retry when the pipeline has finished. | `analysis.jl:617` · §B3 |
| `r.state_leaked` | D | warning | yes / no | no | `R objects from a failed fit remain in the global environment for the life of the process.` | None at runtime; fixed by cleaning up in `finally`. | `estimation.jl:478` · §B2 |
| `r.warning_uncollected` | T | notice | yes / no | no | `R emitted {n} warning(s) outside the per-feature handler; they are in the process log only.` | Grep the log by `run_id`. | `estimation.jl:741` · §B5 |
| `r.adapter_packages_unpinned` | D | warning | yes / no | no | `The adapter records r_packages {packages}, which renv.lock does not pin; nothing verified them.` | Record what is actually loaded. | `Execution.jl:133` · §B4 |

## `ord.*` / `analysis.*` — ordination, diversity, charts

| id | cat | sev | usable / blocks | banner | summary | action | site · audit |
|---|---|---|---|---|---|---|---|
| `ord.nmds.not_run` | R | danger | **no / yes** | yes | `NMDS was not run: {reason}.` | Fix the cause; the coordinates are absent, not zero. | `analysis.jl:1020` · §B3 |
| `ord.nmds.nan_fallback` | R | danger | **no / yes** | yes | `run_nmds returned a full-size matrix of NaN rather than a not-run status.` | None at runtime; fixed by a typed outcome. | `analysis.jl:1020` · §B3 |
| `ord.permanova.not_run` | R | danger | **no / yes** | yes | `PERMANOVA was not run: {reason}.` | Fix the cause. | `analysis.jl:1086` · §B3 |
| `ord.permanova.no_covariates` | Q | notice | yes / no | no | `No metadata column has two or more distinct values, so there is nothing to permute against.` | Supply a grouping. | `analysis.jl:1057` · §B3 |
| `analysis.samples_dropped_below_depth` | Q | notice | yes / no | no | `{dropped} of {total} samples fell below the resolved depth {depth} and were excluded.` | Note the reduced n. | `diversity.jl:131` |
| `analysis.significance.not_run` | R | danger | no / yes | yes | `{test} was not run: {reason}.` | Fix the cause; the chart is unannotated, not annotated with a null result. | `analysis.jl:605-622` |

## `bench.*` — benchmarks

| id | cat | sev | usable / blocks | banner | summary | action | site · audit |
|---|---|---|---|---|---|---|---|
| `bench.allocation_churn` | T | notice | yes / no | no | `{workload} at {taxa} taxa allocated {bytes} (churn, not footprint) in phase {phase}.` | Compare churn against the phase baseline. | `bench/ilr_bases/benchmark.jl:118` · §W4 |
| `bench.peak_rss_growth` | T | notice | yes / no | no | `{workload} grew peak RSS by {bytes}; note maxrss never decreases, so a reusing workload reads 0.` | Read the two metrics together. | `bench/ilr_bases/benchmark.jl:119` · §W4 |
| `bench.duration` | T | notice | yes / no | no | `{workload} took {seconds}s.` | Compare on the same runner only. | `bench/ilr_bases/benchmark.jl:117` |

## `ci.*` / `provenance.*` / `pipeline.*`

| id | cat | sev | usable / blocks | banner | summary | action | site · audit |
|---|---|---|---|---|---|---|---|
| `ci.no_verdict` | D | fatal | **no / yes** | yes | `CI produced no verdict for this commit ({n} consecutive startup failures, zero jobs).` | Fix the workflow; treat the code as untested until then. | `.github/workflows/ci.yml` · §B1 |
| `provenance.component_unproved` | D | warning | yes / no | yes | `{component} could not be proved; this run is degraded: {reason}.` | Restore the component. | `provenance.jl:742` |
| `provenance.database_unproved` | D | warning | yes / no | yes | `Database '{name}' could not be proved; this run is degraded: {reason}.` | Restore or re-pin the database. | `provenance.jl:755` |
| `provenance.hash_sidecar_unreadable` | T | info | yes / no | no | `A hash sidecar could not be read and was ignored.` | None. | `provenance.jl:118` |
| `pipeline.taxa_column_remapping_skipped` | Q | warning | yes / no | no | `Column remapping was skipped: {columns} are not in the data.` | Check the taxonomy table. | `merge_taxa.jl:322` |

---

## Coverage of the audit's findings

Every finding in the audit maps to at least one id, or to a code change with no
user-visible notice:

| audit finding | notice id(s) |
|---|---|
| M1 quadratic default Helmert | `exec.memory.allocation_churn` |
| M2 dendrogram peak is 2× | `ilr.dendrogram.variation_matrix_size` (wording corrected) |
| M3 `clr_table`, M4 copies, M6 DuckDB hand-off, M7, M8 | no notice — code changes only |
| M5 SBP materialisation | `ilr.source.materialised` |
| N1 mixed healing scales | `exec.heal.nan_inf.value_scale_mixed` |
| N2 healing changes data | `exec.heal.*`, `exec.heal.positions_unrecorded` |
| N3 stub check | `exec.diag.batch_confounding_not_implemented` |
| N4 two `rarefy`s | `exec.rarefy_is_scaling` |
| B1 CI | `ci.no_verdict` |
| B2 R state leak | `r.state_leaked` |
| B3 not-run inconsistency | `ord.*.not_run`, `ord.nmds.nan_fallback`, `analysis.significance.not_run` |
| B4 `MASS` / adapter packages | `r.package_not_probed`, `r.package_missing`, `r.adapter_packages_unpinned` |
| B5 uncollected R warnings | `r.warning_uncollected` |
| B6 inert `JuliaAdapter` | no notice — documentation |
| W1 `:not_run` reported safe | `exec.estimation_not_run` |
| W2 disabled analyses unrecorded | the whole `r.*` / `ord.*` set |
| W3 free text | this catalogue |
| W4 benchmark conflation | `bench.allocation_churn`, `bench.peak_rss_growth` |
| W6 repeated warnings | `r.unavailable` (memoised) |
