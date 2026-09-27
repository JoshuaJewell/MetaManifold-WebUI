<!--
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Baton handoff — MetaManifold-WebUI issue #1, formal verification

Paste everything below the line to the next agent.

---

You are taking over work on `hyperpolymath/MetaManifold-WebUI`, GitHub issue #1
("Add a validated Julia statistics layer with exact counts/rationals and higher
precision"). The Agda formal-verification layer is **built, verified, pushed, and
open as draft PR #79**. Your job is to land it and then continue the validation
work. Read this whole document before touching anything — the first three
sections contain traps that already cost a previous session significant time.

## 1. The state of the repository (verify this before trusting it)

**`main` has been replaced and has NO common ancestor with the working branch.**

```sh
git fetch origin main
git rev-list --count origin/main      # -> 1
git merge-base HEAD origin/main       # -> empty
```

`origin/main` is `770a614` ("feat(proofs): Evidence Mode foundations — Agda
library + golden vectors (#7) (#78)"), a single-commit history. Any notes,
session memory or plan that says main is at `e33e10f` are stale.

`main` **already contains an Agda suite**, on a different toolchain:

```
proofs/agda/MetaManifold/All.agda
proofs/agda/MetaManifold/Composition/{Tree,Node,Sum}.agda
proofs/agda/MetaManifold/ILR/{SBP,Contrast,Kernel,Invariance,Orthonormal,Comb,Integer}.agda
proofs/agda/MetaManifold/Evidence/{Prelude,Residual,Echo,Warrant,Signed,Decision}.agda
proofs/agda/reject/{Postulate,OrthogonalSelf,KernelWithoutHypothesis,IdentificationWithoutUniqueness}.agda
proofs/agda/metamanifold-proofs.agda-lib
proofs/vectors/evidence_vectors.json
docs/formal/verification-plan.md
docs/formal/evidence-verification.md
```

Read `docs/formal/verification-plan.md` on `main` **first**. It defines the
estate's scope and residue conventions and your work has to fit them.

## 2. Three traps that already bit

1. **`main` pins Agda 2.6.4.3 / agda-stdlib 2.1. The new proofs do not compile
   there.** They were written against agda-stdlib's development line towards 3.0
   (`name: standard-library-3.0`), pinned by SHA
   `2ffa8b7d4e8e818717ad643d184f055a4d1b0447`. **No `v3.0` tag exists** — the
   newest is `v2.4` — and against `v2.4` the first failure is
   `Data.Integer.Properties` not exporting `_≡?_`. Porting targets, in the order
   they break, are in `proofs/residue/toolchain.residue` (`R-TC-1`). Note that
   `main`'s `Evidence.*` modules are deliberately **stdlib-free**
   (`Agda.Builtin.*` only) so they check under any Agda ≥ 2.6.4.3 — matching that
   design removes the porting problem entirely, since `ℚᵘ` is the only thing in
   the statistics modules that genuinely needs the stdlib.

2. **Four concrete collisions block a naive merge.** Both sides define
   `module MetaManifold.All`; both have a `.agda-lib` with `include: .` in the
   same directory (mine `metamanifold`, main's `metamanifold-proofs`); `main`'s
   `Justfile` already defines a `proofs:` recipe at line 279; and the toolchains
   differ. Full detail is in the first comment on PR #79.

3. **`/tmp` does not survive between turns in the Arena sandbox.** A toolchain
   installed there silently vanishes and the proofs then cannot be re-verified.
   `proofs/bootstrap.sh` installs into `proofs/.vendor/` (git-ignored) for this
   reason. Do not "fix" it back to `/tmp`.

## 3. What exists and is verified

Branch `arena/01a0db37-metamanifold-webui`, HEAD `b9694a0`, PR #79 (**draft**).

Six modules, 77 top-level definitions, plus a gate entry `MetaManifold/All.agda`:

| Module | The substance |
| --- | --- |
| `Prelude` | `Outcome A = value A ⊎ refused Refusal`, `outcome-total` (no third arm), `ℚᵘ` helpers, `whole-is-one` with **no `n ≢ 0` hypothesis** because `ℚᵘ`'s denominator is `suc _` — a zero denominator is unrepresentable by type. This is why Agda and not Lean. |
| `Proportions` | `relativeAbundance` refuses in the order the Julia layer checks; proportions sum to exactly `1ℚᵘ`; a returned value *is* `c/t`. |
| `ExactCounts` | `checkedAdd` exact **in both directions** (returned value is the true sum; refusal happens exactly when the sum leaves the range) + `checkedSumOf-is-exact` for the fold. |
| `PermutationTest` | `never-reports-zero` for every input, plus a proved negative control showing the replaced estimator really does return exactly zero. |
| `BenjaminiHochberg` | scaling is exactly `(M/(j+1))·(n/(d+1))`; a bigger family never makes a q-value smaller. |
| `DecimalRounding` | "correctly rounded to *s* decimals" as a predicate over integers alone; machine-checked known-answer vectors; tie handling stated in both directions and proved to agree elsewhere. |

Gate machinery, all verified **from an empty `proofs/.vendor`** (the bootstrap
installed its own toolchain):

- `proofs/bootstrap.sh` → `proofs: OK`, exit 0
- `proofs/tests/axiom-audit.sh` → `7/7 modules reachable`, `clean`, exit 0
- `proofs/tests/gate-selftest.sh` → `10/10 controls behaved correctly`, exit 0
- `.github/workflows/proofs.yml` → **has never executed.** Neither has `ci.yml`
  on this PR, and `gh api .../actions/workflows` does not list `proofs.yml` at
  all, so the cause is repo/org-level rather than the file. Check this early.

`proofs/tests/gate-selftest.sh` is the piece with no counterpart on `main`: it
breaks the proofs nine ways on purpose (wrong known-answer digit, reversed
monotonicity, plus-one removed, zero total silently zero, overflow no longer
refused, module dropped from the gate entry, postulate injected, `--safe`
removed, tie forced the wrong way) and requires each to be rejected. `main`'s
`proofs/agda/reject/` expresses the same idea differently — they are
complementary, not duplicative. On its first run the self-test reported **1/10**
because `bash -c "$mutator" "$file"` binds the path to `$0`; the control reported
"mutation did not apply" instead of passing. That is the point of having it.

## 4. What is NOT proved

`proofs/residue/` has every open obligation with an id, a precise statement, an
impact rating and what would close it. Read it before claiming anything.
Highlights:

- `R-EC-1` `Checked.refused` injectivity — open, low impact.
- `R-BH-1` non-negativity of the scaled value — open, mechanical obstacle
  recorded (`+ M * + n` does not reduce to `+ (M * n)`).
- `R-BH-2` the step-down `min` envelope — not attempted; `ℚᵘ` has no `⊓` lemmas.
  **`q_i ≥ p_i` is FALSE in general — do not attempt it.**
- `R-TC-1` the toolchain port. `R-TC-2` (closed) Agda's library-file location.
- `out-of-scope.residue` — **no probability theory, no IEEE-754, no proof that
  `numeric_policy.jl` implements the model, nothing about the `:ordinary` path.**

## 5. Your task, in order

1. `git fetch origin main` and read `docs/formal/verification-plan.md` and
   `proofs/agda/README.md` on `main`. Confirm §1 above still holds.
2. Rebase `arena/01a0db37-metamanifold-webui` onto `origin/main`:
   - move the six modules to `MetaManifold/Statistics/*`
   - delete `metamanifold.agda-lib`; add the modules to main's `All.agda`
   - rename the `Justfile` recipes (`proofs-statistics-*`) or fold them into
     main's existing `proofs:` recipe and `scripts/check-proofs.sh`
   - **either** port to stdlib 2.1 **or** make the modules stdlib-free like
     `Evidence.*` (preferred; see `R-TC-1`)
3. Re-run the gate and **re-verify**. Do not assume a port is correct because it
   was correct before.
4. Investigate why no in-repo workflow ran on PR #79. Un-draft only when the
   proofs workflow has actually been observed to pass on a runner.
5. Then continue issue #1 itself — the Julia validation layer. **Julia is not
   installable in the Arena sandbox** (no julialang-s3, no
   `objects.githubusercontent.com`, so `juliaup`/`elan` are dead; no
   `pkg.julialang.org`; `apt-get` unusable). Anything Julia you write there must
   be explicitly marked unrun and wired into CI. Still open on #1: per-method
   known-answer + independent-reference validation, simulation studies with
   pre-set tolerances, adversarial/property + fuzz coverage,
   UI-to-backend-to-result workflow/accessibility/provenance tests, benchmarks,
   existing-behaviour-unchanged verification, concurrency/cancellation/
   resource-limit tests, CI standards gates with retained artifacts, and
   independent statistical review.
6. `test/fixtures/agda-known-answers.json` is the Agda→Julia bridge: vectors the
   proof assistant checked against a definition, each annotated with the lemma
   that fixes it. Write the Julia conformance testset against it. If the two
   disagree, the Julia layer is wrong.

## 6. Agda syntax traps (all confirmed the hard way)

- Never write an operator with underscores in expression position (`a _*_ b`) —
  it parses as a name/section and the error names only the *other* operators.
- Fixity does not travel with `renaming`. Import unqualified, rename the other
  operator.
- Once ℤ's `_≤_` and any rational `_≤_` are both in scope, every bare `≤` is
  ambiguous. Qualify one.
- `open IsEquivalence … renaming …` / `open <Solver> …` do not re-export — add
  `public` in the Prelude.
- `open import` never re-exports, so each module imports `_≃_`, `mkℚᵘ`, `yes`/`no`,
  `⊥-elim` itself.
- Nested `with` is fragile; split into separate top-level clauses.
- `λ ()` only discharges a goal of type `⊥`.
- `List.map f (List.map g xs)` is not definitionally `List.map (f ∘ g) xs`.
- `n * 1` and `1 * n` on ℕ do not normalise.
- Proving properties of a `with`-based function is far easier if the decision is
  returned as **data** (`AddOutcome` in `ExactCounts`) and matched on.
- `--without-K` rejects `rewrite` when the resulting index is not a variable.
- `solve n f refl` is ∀-quantified — apply the variables *after* `refl`.

## 7. Environment facts

- Verified invocation: `proofs/bootstrap.sh` (or
  `proofs/bootstrap.sh --check` once bootstrapped), from `proofs/agda`.
- Usable network from the sandbox: GitHub code/clone/API (when the token is
  valid), PyPI, npm. **Not** usable: julialang-s3,
  `objects.githubusercontent.com`, `pkg.julialang.org`, `deb.debian.org`,
  hackage/ghcup, Lean hosts.
- `pip3 install` on the system interpreter fails (PEP 668); `apt-get` is unusable.
- `gh issue view <n>` without `--json` fails (Projects classic deprecation).
- GitHub auth in the sandbox **expires between turns**. If `git push` fails with
  `could not read Username`, ask the user to reconnect GitHub in Arena.
- Commits can be lost between turns if HEAD is reset; always re-check
  `git rev-parse HEAD` before claiming something is pushed, and confirm with
  `git ls-remote origin refs/heads/<branch>`.
