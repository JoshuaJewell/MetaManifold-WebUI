<!--
SPDX-License-Identifier: CC-BY-SA-4.0
-->
# Integrating the fork into upstream — a low-friction path

**Audience:** Joshua (upstream maintainer) and the `hyperpolymath` engineering
team. **Goal:** get the fork's work into
`JoshuaJewell/MetaManifold-WebUI` without a single "hundreds of files, load of
conflicts" review, and without either side ever being unable to run.

The measured problem and the exact file counts are in
[`conflict-map-2026-09-25.md`](./conflict-map-2026-09-25.md). In one sentence:
**the fork and upstream share only the root commit, diverged at the second
commit, and then developed in parallel for months — so the merge base is nearly
empty and a normal merge surfaces all 159 shared‑but‑different files at once.**
Nothing in this guide changes upstream's operational behaviour unless you choose
it to; the default is "behave exactly like upstream".

For the separate KYAML authoring-format pilot, see the [upstream architecture and
transition audit](./kyaml-upstream-adoption.md). It identifies the runtime config readers,
user-editable YAML writers, GitHub Actions boundary, and the evidence required before
changing any tracked files.

## The idea, in layers

0. **Re‑anchor** — replay the fork's commits one‑by‑one onto upstream so git
   auto‑applies the clean ones and only the genuine overlaps surface, per commit.
   Measured: the whole fork re‑anchors with **3 decisions** and **0 residual
   conflicts**, producing ~206 granular, reviewable commits. (`just reanchor`.)
1. **Merge hygiene** — stop git from ever producing a line‑conflict inside a
   lockfile or a build artefact. (`.gitattributes` + `just merge-drivers`.)
2. **Profiles & components** — turn "review 159 files" into "make ~7 ordered
   trust decisions", switchable from *pure upstream* to *partially transitional*
   to *everything*. (`config/integration.toml` + `scripts/integrate.sh` + the
   `just integrate …` recipes.)
3. **A staging plan** — land the safe infrastructure first, let CI go green, then
   take the risk‑bearing pieces one tier at a time. (`just integrate-plan`.)

Everything below is additive and default‑off. A fresh clone behaves as
`base` — i.e. exactly like upstream — until someone runs a command to change it.

---

## 0. Re‑anchor — from one 159‑file wall to ~206 granular commits

This is the deepest fix and the reason the fork *looks* unmergeable. Instead of
one giant three‑way merge, re‑anchor replays the fork's history onto upstream:

```bash
just reanchor-plan      # read-only: classify every fork commit (auto / overlap / delete-risk)
just reanchor           # do it; leaves branch reanchor/onto-upstream for review
just reanchor-manual    # same, but stop at every conflict for hands-on resolution
```

It runs in an isolated worktree and never touches your current branch. The result
is upstream's history with the fork's ~206 commits cleanly on top — reviewable and
mergeable incrementally, and sharing a real base so future upstream work merges
cleanly. See `scripts/reanchor.sh` and the measured numbers in
[`conflict-map-2026-09-25.md`](./conflict-map-2026-09-25.md).

---

## 1. Merge hygiene — never hand-merge a generated file

The committed `.gitattributes` now marks lockfiles and build output so git will
not line-merge them (`merge: unset`, `linguist-generated`). Concretely, these can
never again appear as a wall of `<<<<<<<` markers:

`Manifest.toml`, `renv.lock`, `frontend/bun.lock`, `bun.lockb`,
`package-lock.json`, `pnpm-lock.yaml`, `renv/activate.R`, `renv/library/**`,
`config/tools.yml`, `web/dist/**`, `frontend/dist/**`, coverage, `*.cov`,
`lcov.info`, test-results, `*junit.xml`, `bench/**/results/*.json`.

These are resolved by **regenerating**, never by editing:

```bash
just heal        # re-syncs pins, re-instantiates Julia, restores renv, reinstalls frontend
```

If you want lockfiles to *auto*-resolve to the branch you are merging into (and
then be regenerated), wire the stronger opt-in driver once:

```bash
just merge-drivers   # writes merge.lockfile into local .git/info/attributes (never committed)
```

This is fully opt-in and local: it cannot break a merge on a clone that has not
run it, and it never mutates a tracked file.

## 2. Profiles & components — the transitional dial

`config/integration.toml` is the shared registry. It defines the fork's
contributions as named **components**, each with an area, a risk level, the path
globs that identify its files (used to auto-classify conflicts), the upstream PRs
that carry it, and — where the code supports it — the runtime gate that keeps it
safe while transitional. It defines three **profiles** that bundle components:

| profile | what is live | use it when |
|---|---|---|
| `base` | nothing risk-bearing from the fork | you want to behave exactly like upstream (the default) |
| `transitional` | safe infra only: estate tooling, download retries, the 1×1 fix, benchmarks/CI, exact offsets | you want the reliability wins but not the statistics/UI changes yet |
| `full` | everything, including the real-statistics engine and the UI migration | you have reviewed and trust the higher-risk pieces |

Inspect and switch:

```bash
just integrate status                 # active profile + each component on/off
just integrate profiles               # what each profile bundles
just integrate profile transitional   # switch (machine-local, never committed)
just augment  real-statistics         # turn one component on for this checkout
just suspend ui-migration             # turn one component off
```

The active choice is **machine-local** (`METAMANIFOLD_INTEGRATION_PROFILE`, or
the git-excluded `config/integrate.active`). Your choices never become a merge
conflict of their own.

### Why this never compromises operations

- The default profile (`base`) activates nothing; behaviour is unchanged.
- Every component is either a **pure addition** (tooling, CI, retries) or is
  **gated by the code's existing refusal architecture**. The statistics engine
  already refuses loudly rather than guess: a method is either wired to a real
  model (`SUPPORTED_DISPERSION`) or it refuses (`REFUSED_DISPERSION`) with an
  explicit, actionable error. That refusal *is* the transitional state — a method
  can stay refused until it is trusted, so there is no path to a silently-wrong
  result. See `docs/statistics/` and `docs/owner-review-2026-09-25.md`.
- `just heal` and `just doctor` only repair or report; they change no scientific
  behaviour.

## 3. The staging plan

```bash
just integrate plan
```

prints the recommended order (safest → riskiest), cross-referenced with the
upstream PR numbers:

1. `estate-tooling` (none) — Justfile, mise, `.gitattributes`, hooks, check
   scripts. Pure additions; upstream has none of these, so **zero conflict**.
2. `install-retry` (low, PR #11) — retry/TLS around byte-exact downloads.
3. `matrix-1x1-fix` (low, PR #8) — the 1×1 matrix edge case.
4. `benchmarks-ci` (low, PR #10) — baselines + CI matrix.
5. `exact-offsets` (medium) — TSS/CSS/TMM/RLE offsets + lowercase method match.
6. `real-statistics` (high) — real models + refusals (no invented numbers).
7. `ui-migration` (high, PR #7) — the reactive-UI modernisation.

Land one tier, let CI go green, then the next. This is upstream's "Option B"
(selective cherry-picking) from the owner review, made mechanical.

---

## Workflow for the maintainer (Joshua)

```bash
# once, in your upstream clone
git remote add fork https://github.com/hyperpolymath/MetaManifold-WebUI.git
git fetch fork

# see exactly what you are taking on, by category
just integrate-triage          # during a merge/rebase; or point it at a list:
just integrate status
just integrate plan

# integrate the safe infrastructure first (tiers 1–4), regenerate locks, gate:
git checkout -b integrate-infra fork/main
#   … resolve only the ~handful of non-generated conflicts the triage flagged …
just heal && just ci

# then opt into the risk-bearing pieces a tier at a time, trusting as you go:
just integrate profile transitional   # try the safe set live
just integrate profile full           # once you trust the statistics + UI
```

During any merge/rebase, `just integrate-triage` classifies every conflicted path
as **AUTO** (regenerate — lockfiles/build output), **COMPONENT** (decide per the
registry: take the side for components you have chosen to trust), or **HUMAN**
(unclassified — review). That collapses "hundreds of files" into "a handful of
human decisions plus a `just heal`".

## Workflow for the fork team — stay turnkey *and* merge-friendly

**Turnkey (your side always runs with zero extra effort):**

```bash
just bootstrap     # mise toolchain (julia/bun/node/just) + frontend deps + hooks + merge drivers + tool map
just setup-full    # …plus Julia instantiate, renv restore, and sha256-pinned external tools
just doctor        # health report (fails loudly if an essential tool is missing)
just heal          # repair anything that drifted
```

`mise.toml` is the single source of truth for the toolchain (every pin verified
against CI); `just setup-tools` provisions it; `just renv-restore` restores the
byte-exact R lockfile; `install.sh` installs the sha256-pinned external pipeline
tools. R is the one documented exception (not in the mise registry) and is called
out by `just doctor` rather than silently skipped.

**Merge-friendly (keep the conflict surface small):**

- Put new work in **new files/directories** wherever possible — a clean addition
  never conflicts.
- Keep edits to shared root files (`start.sh`, `install.sh`, `install.jl`,
  `Project.toml`, `.github/`, `README.md`) minimal and in append-safe regions.
- Never hand-edit lockfiles/build output; run `just heal`/`just build` and commit
  the regenerated artefact.
- Register every contribution in `config/integration.toml` so triage can classify
  it automatically.

## Quick reference

| command | what it does |
|---|---|
| `just doctor` | environment health (hard-fails on missing essentials) |
| `just heal` | repair the environment to a known-good state |
| `just bootstrap` / `just setup-full` | clone-to-runnable on a bare machine |
| `just merge-drivers` | opt-in lockfile auto-resolution (local, safe) |
| `just reanchor-plan` / `just reanchor` | re‑anchor the fork onto upstream as granular commits |
| `just integrate status` / `profiles` / `profile <p>` | the transitional dial |
| `just integrate plan` / `triage` / `verify` | staging order / conflict classes / gates |
| `just augment <id>` / `just suspend <id>` | flip one component |

Engine: [`scripts/integrate.sh`](../../scripts/integrate.sh). Registry:
[`config/integration.toml`](../../config/integration.toml).
