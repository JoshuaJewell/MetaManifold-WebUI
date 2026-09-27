<!--
SPDX-License-Identifier: CC-BY-SA-4.0
SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
-->
# KYAML adoption: upstream architecture and a low-risk path

**Purpose:** explain where an authoring-format change would touch Joshua's upstream
repository, separate safe formatting from runtime behavior, and give each review a small,
explicit scope. This is an architecture review, not a request to adopt KYAML or a claim that
an upstream change has been tested or merged.

**Snapshot reviewed:** `JoshuaJewell/MetaManifold-WebUI` `main` at
`ecefb1c72b3d2515e7086024b14227ef13329605` (2026-09-27). Recheck paths and tests against
the target commit when preparing a PR.

## Bottom line

KYAML is a source-authoring convention, not a replacement parser or a new application data
format. Upstream already reads configuration with `YAML.jl` 0.4. If a file remains valid
YAML and its parsed values are unchanged, its existing consumers need no new dependency and
no API or pipeline change. **Keep `YAML.jl`, existing serializers, and user configuration
formats unchanged during the pilot.**

The boundary that matters is who owns and writes each file:

| Files / surface | Current consumer or writer | Safe transition rule |
|---|---|---|
| `config/defaults/*.yml` | Julia `YAML.load_file`; pipeline defaults feed the default → global → study → group → run cascade in `src/core/config.jl`. | Convert as repository-owned inputs; compare the values read by `YAML.jl` before/after. Do not change the cascade or hashes as part of a format-only patch. |
| `config/ci/*.yml`, `bench/**.yml`, tracked sample `data/MiSeq_SOP/pipeline.yml` | CI/bootstrap, benchmark, or example-data readers. | Convert as authored inputs and run their existing consumers. Keep the sample's values and provenance intact. |
| `.github/workflows/*.yml` | GitHub Actions' workflow parser—not the application YAML loader. Several files contain multiline shell `run:` blocks; Dependabot and `gh actions-lock` also rewrite workflow pins. | Convert only after an actual GitHub Actions parse/run probe. Preserve trigger, job IDs, `needs`, permissions, matrix values, environment, working directory, shell, and `GITHUB_ENV` behavior. Bot drift may exempt these files from the canonicality gate, but not from conversion or rollback. |
| `config/pipeline.yml`, `config/primers.yml`, `config/databases.yml`, `config/composition.yml`, `config/tools.yml`, and per-study `pipeline.yml` overrides | Machine/user-editable overlays; routes in `src/server/routes/config.jl`, `composition.jl`, and `results.jl` write YAML with `YAML.write`. Root `.gitignore` excludes the machine-level config files. | Do not reformat user state or change save behavior in the pilot. The reader continues to accept YAML; runtime writes may remain block-style YAML. Keep generated/user files outside the tracked-file gate. |
| `run_config.yml`, provenance sidecars, result/filter files | Generated from the config cascade or written as run output; `src/core/config.jl` and `src/core/provenance.jl` own these paths. | Leave generated artifacts and their writers alone. They are not repository-authored source configuration. |

This means Joshua does **not** have to take a KYAML dependency or immediately migrate
user-generated files. A format-only change to the tracked defaults can be read by the
existing YAML loader; an edited user overlay can continue to be written in ordinary YAML.

## Risks that deserve explicit evidence

1. **Semantic equivalence, not visual similarity.** Quoting is intentional: KYAML quotes
   strings and keys that YAML readers may otherwise interpret differently (`on`, `no`,
   version-like strings). The gate must compare the actual parsed data from the production
   `YAML.jl` reader before and after emission for every tracked configuration document.
   The emitter's own YAML→KYAML→YAML round-trip is not enough by itself.
2. **GitHub workflow interpretation.** The quoted `"on"` key must still trigger the
   workflow. Unit tests cannot establish this. A real PR run must report the expected jobs
   under their unchanged stable names, with no missing required status context.
3. **Shell extraction is a separate behavior change.** Moving `run: |` blocks to
   `scripts/ci/*.sh` can change the working directory, shell options, expression expansion,
   environment, exit behavior, permissions, or writes to `GITHUB_ENV`. Keep that work in a
   separate reviewed change before the workflow-format rewrite, with a shell-level test or
   the same real CI execution.
4. **Bot rewrites.** A bot reverting a workflow to block style should not make the quality
   gate flaky, but operators must still be able to canonicalize and roll back those files.
   The exemption must apply only to the gate—not to the converter.
5. **User data and serializer scope.** Applying the converter to ignored machine configs,
   run directories, or provenance outputs would create churn and could rewrite user state.
   The source-format pilot should enumerate tracked YAML only and must not alter the
   application's YAML serialization paths.
6. **Reversibility.** Preserve a byte-exact Git revert for the format-only conversion.
   `just use-yaml` is a canonical re-emission, not a promise to recover every original
   whitespace byte; the Git revert is the exact escape hatch.

## Small, reviewable sequence

### 0. Hardening patch — no source YAML rewritten

Land the KYAML tool correction and its regression tests first. The drift list is a gate
allowlist, not a migration exclusion. Test the converter and rollback against an exempt
workflow path, and parse every tracked YAML document. Compare application-consumed
configuration with the production `YAML.jl` reader; validate workflow files with GitHub
Actions itself. Correct the corpus census and document the actual mutable-file boundary.
This patch changes no workflow or application behavior.

### 1. Fix and establish the CI verdict

Resolve the current fork CI `startup_failure` before treating CI as evidence. A run with
zero jobs has no test verdict. Then establish a passing baseline for the KYAML unit tests,
workflow syntax checks, and existing full test suite on the exact review commit.

### 2. Extract CI shell blocks (only if still required)

If the project retains the pilot's script-extraction decision, do it separately from
formatting. For every moved step, record the old and new `working-directory`, shell,
`env`, GitHub expression, permissions, `set -e`/pipefail, and outputs. Run the workflow
before changing its YAML style.

### 3. Format-only conversion

Convert only tracked, repository-authored `.yml`/`.yaml` files. Keep the migration as one
format-only change (or one clearly named commit within the PR) so its purpose and revert
boundary stay obvious. Its review description should report:

- the exact file list and parser-equivalence result;
- every decision count from `just kyaml-report` (`null`, quoting, comment movement);
- the GitHub Actions run IDs that parsed and executed the converted workflows;
- the required status contexts observed;
- files intentionally not gated because of bot rewriting; and
- confirmation that `Project.toml`, serializers, APIs, config-cascade behavior, run hashes,
  and user data are unchanged.

Enable the canonicality gate in the same migration change. Do not claim the GitHub workflow
proof is complete merely because the ordinary YAML parser accepts the file.

## Findings from this checkout that must be resolved before migration

- The repository has **17 tracked YAML files**, not 16. Two workflow paths are listed in
  `config/kyaml/drift.txt`; that makes 15 files subject to the canonicality gate after
  conversion, while all 17 still need parser/round-trip coverage.
- Before the current hardening change, the CLI filtered the drift paths out of every mode.
  That meant the workflows would not be converted by `use-kyaml` and would not be restored
  by `use-yaml`, despite the pilot saying they are in scope. The implementation now applies
  those exemptions only to `--check`; a regression test covers conversion, gate exemption,
  and rollback.
- The CLI also referred to `Stats` from outside its `KYAML` module without qualifying the
  name. It now uses `KYAML.Stats`; the new CLI-level regression test exercises that path.
- The old parser-corpus test skipped the two workflows. It now includes every tracked file,
  compares both emitted forms with `YAML.jl` for application-consumed YAML, and leaves
  workflow acceptance to GitHub Actions itself. The current Actions workflow contains only
  a placeholder comment about adding the canonicality gate; it does not yet run
  `check-kyaml`. Add the real gate in the same change that converts the corpus.
- **These Julia tests have not yet been executed in this sandbox**: Julia is absent, and the
  latest fork CI run ended in `startup_failure` before starting any jobs. Therefore this
  audit does not certify that all 17 files pass or that GitHub has accepted a converted
  workflow.

## Upstream source pointers

- [`Project.toml`](https://github.com/JoshuaJewell/MetaManifold-WebUI/blob/main/Project.toml) — YAML.jl 0.4 remains a runtime dependency.
- [`src/core/config.jl`](https://github.com/JoshuaJewell/MetaManifold-WebUI/blob/main/src/core/config.jl) — defaults and override cascade; generated `run_config.yml`.
- [`src/server/routes/config.jl`](https://github.com/JoshuaJewell/MetaManifold-WebUI/blob/main/src/server/routes/config.jl) — config reads and atomic YAML writes for editable overlays.
- [`src/core/primers_library.jl`](https://github.com/JoshuaJewell/MetaManifold-WebUI/blob/main/src/core/primers_library.jl) and [`src/core/databases_library.jl`](https://github.com/JoshuaJewell/MetaManifold-WebUI/blob/main/src/core/databases_library.jl) — user-edited library files and normalization/validation boundaries.
- [`.github/workflows/ci.yml`](https://github.com/JoshuaJewell/MetaManifold-WebUI/blob/main/.github/workflows/ci.yml) — CI parser, shell steps, Julia matrix, and tool installation.
