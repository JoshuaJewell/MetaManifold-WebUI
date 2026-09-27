<!--
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Proof status

Formal verification of the validated Julia statistics layer (issue #1), in **Agda
2.7.0.1** with **agda-stdlib pinned at `2ffa8b7d4e8e818717ad643d184f055a4d1b0447`**,
under `--safe --without-K`.

> **Read this before trusting the pin.** That SHA is on agda-stdlib's development
> line towards 3.0 — the branch's own library file declares
> `name: standard-library-3.0`, but **no `v3.0` tag exists** (the newest is
> `v2.4`). These proofs do **not** compile against `v2.4`, nor against `v2.1`,
> which is the estate pin on `main`. See `R-TC-1` in
> [`residue/toolchain.residue`](residue/toolchain.residue).

| Gate | Command | Result |
| --- | --- | --- |
| Type-check every proof | `proofs/bootstrap.sh --check` | **PASS** — exit 0, no warnings |
| Axiom audit | `proofs/tests/axiom-audit.sh` | **PASS** — 7/7 modules reachable, no postulates, no FFI, no unsound flags, no holes, all `--safe` |
| Gate self-test (does the gate reject anything?) | `proofs/tests/gate-selftest.sh` | **PASS** — 10/10 controls |

Last verified locally: 2026-09-26, from an **empty `proofs/.vendor`** — i.e. the
whole toolchain was installed by `proofs/bootstrap.sh` itself (PyPI wheel for
Agda, stdlib fetched at the pinned SHA, Agda source for the primitive libraries)
and then:

- `proofs/bootstrap.sh` → `proofs: OK`, exit 0
- `proofs/tests/axiom-audit.sh` → `7/7 modules reachable`, `clean`, exit 0
- `proofs/tests/gate-selftest.sh` → `10/10 controls behaved correctly`, exit 0

## Why Agda and not Lean

The brief allowed Lean as a fallback. Agda was used because it is the right tool
for *this* subject rather than the merely available one: the claims are about
exact rational arithmetic and about total functions that either return a value or
name a refusal. Dependent indices let `ℚᵘ`'s denominator be `suc denominator-1`,
which makes **a zero denominator unrepresentable by type** — the refusal for
division by zero is not a checked branch, it is an impossibility. That is a
stronger statement than a guarded division and it comes for free.

## Modules

Every module is reachable from `MetaManifold/All.agda`; the axiom audit fails the
build if one is not.

### `MetaManifold.Prelude` — 22 definitions
The shared model: `Refusal` (five constructors, each annotated with the Julia
exception it stands for), `Outcome A = value A ⊎ refused Refusal`, exact rational
helpers over `ℚᵘ`, and the summation lemmas.

Key results:

- `outcome-total` — every outcome is either a value or a named refusal. There is
  no third arm; this is a statement about the shape of the type, so a later edit
  that forgets a case cannot violate it quietly.
- `refused≢value`, `value-injective`, `is-value`
- `whole-is-one : ∀ n → mkℚᵘ (+[1+ n ]) n ≃ 1ℚᵘ` — **with no `n ≢ 0`
  hypothesis**, because `ℚᵘ` cannot represent a zero denominator. This is the
  exact sense in which the refusals are total.
- `zero-over`, `over-denominator`, `common-denominator-+`
- `sumℕ`, `sumℤ`, `sumℚᵘ`, `sumVecℕ`, `sumVecℤ`, `sumVecℚᵘ`, `sumVecℕ-map-+`
- `sum-over-common-denominator`, `sumVec-over-common-denominator`
- `numerator-monotone-≤` — the single route for comparing two rationals over a
  common denominator. Hand-built cross-multiplication is banned by convention
  because it is where sign errors live.
- `_≟ℚᵘ_`, `_/[_]_`, `+-cong-≃ʳ`, `≡⇒≃`

### `MetaManifold.Proportions` — 13 definitions
Models `exact_relative_abundance` and the exact summary path.

- `relativeAbundance c t : Outcome ℚᵘ` — refusals are checked **in the order the
  Julia layer checks them**, so the reason a caller sees here is the reason the
  Julia layer would give.
- `zero-total-is-refused` — `relativeAbundance 0 zero ≡ refused zeroTotal`.
- `zero-total-with-count-is-refused-as-exceeding` — `relativeAbundance (suc c)
  zero` refuses as `countExceedsTotal`, **not** `zeroTotal`. This was a finding,
  not an assumption: `exact_relative_abundance` checks `count <= total` before
  `iszero(total)`, so a non-zero count against a zero total reports the count
  error first. A test written from the docstring would have expected the other.
- `count-exceeds-total-is-refused`
- `value-is-the-quotient` — when a value comes back, it *is* `c/t` exactly.
- `value-has-nonzero-denominator`
- `abundances`, `abundances-list`, `sum-abundances-list`
- `proportions-sum-to-one`, `proportions-sum-to-one-list` — exact proportions of a
  sample sum to `1ℚᵘ`, with no rounding anywhere in the chain.
- `aggregate-then-divide`

### `MetaManifold.ExactCounts` — 12 definitions
Models `require_count` / `checked_count_sum`: a fixed-width accumulator whose
bound is *data* (`Fits k v = -2^k ≤ v < 2^k`, `pow2`), not an assumption.

- `fits?` — the guard is decidable.
- `checkedAdd`, `checkedSumOf` — the true sum if it still fits, `refused
  overflow` if it does not. There is no branch in which a wrapped value escapes.
- `checkedAdd-is-exact` — when it returns a value, that value is the true sum. It
  never returns a plausible-looking wrong total.
- `checkedAdd-never-refuses-a-sum-that-fits` — when the true sum fits, it returns
  it. The guard is not so cautious that it refuses valid work.
- `checkedAdd-refuses-exactly-when-it-must` — it refuses precisely when the sum
  leaves the range, so nothing wraps.
- `checkedSumOf-is-exact` — the fold version: if a table's counts sum without
  overflowing, the reported total *is* the sum of the counts, and every exact
  proportion downstream inherits that.

Failing in either of the first two directions would be a defect; only the pair
together is the property the layer claims, which is why both are proved rather
than one of them being called "the safety property".

### `MetaManifold.PermutationTest` — 11 definitions
- `pValue b B = (b+1)/(B+1)` — the plus-one (Phipson & Smyth /
  Davison–Hinkley) estimator.
- `naivePValue b B = b/B` — the estimator it replaces, kept on purpose: a
  contrast against nothing proves nothing. Its denominator is `B`, so `B = 0` is
  a refusal, which is itself part of the contrast.
- `never-reports-zero` — the reported p is strictly positive for every input.
  This is what makes "p = 0" unprintable *by construction* rather than by
  convention: no code path, no rounding rule and no display layer has to be
  trusted to avoid it.
- `at-most-one` — it never exceeds 1, so it cannot be read as a likelihood ratio.
- `resolution-is-one-over-B+1` — `1/(B+1)` *is* the Monte Carlo resolution.
- `resolution-is-a-lower-bound`, `monotone-in-extremes`
- `naive-estimator-can-report-zero` — the negative control: the replaced
  estimator really does return exactly zero. A proof that a hazard was avoided is
  only evidence if the hazard is shown to be real.
- `plus-one-does-not-at-the-same-input`, `exhaustive-is-exact`

### `MetaManifold.BenjaminiHochberg` — 4 definitions
- `bhScale M j n d` — `M/(j+1) · n/(d+1)` built by a single `mkℚᵘ`, with the
  denominator expanded, so there is no intermediate rounding.
- `bhScale-numerator`, `bhScale-denominator` — the numerator of the result is
  `M·n` and its denominator is `suc (j + d + j·d)`, i.e. `(j+1)(d+1)`. The `suc`
  is the point: it is where an off-by-one would hide.
- `monotone-in-family-size` — a bigger family can never produce a smaller
  q-value. A multiplicity correction that went the other way would reward running
  more tests, and would do so silently. Takes non-negativity of the numerator as
  a hypothesis, which is not a formality: the scaling multiplies the numerator,
  so a negative numerator would reverse the inequality.

### `MetaManifold.DecimalRounding` — 15 definitions
Specifies what "correctly rounded to *s* decimal places" means, as a predicate
over integers only — no division, so no rounding inside the specification of
rounding.

- `IsRounding n d s m` ≡ `2n·10^s − (d+1) ≤ 2m(d+1) < 2n·10^s + (d+1)`
  (exact halves **down**: `1/2` at 0 dp is `0`).
- `IsRoundHalfUp n d s m` ≡ `2m(d+1) − (d+1) ≤ 2n·10^s < 2m(d+1) + (d+1)`
  (exact halves **up**: `1/2` at 0 dp is `1`).

Which one the display layer must use is a contract decision; that the two differ,
and exactly where, is what these definitions record. Nothing about tie handling is
hidden in a library default.

- `isRounding?` — the spec is decidable, so a harness can ask it of any candidate
  instead of trusting a second implementation of the same rule.
- Machine-checked known-answer vectors (each closed by `toWitness` applied to a
  decision procedure — if the digits did not satisfy the spec, the definition
  would not type-check): `rounding-2/3-at-2dp` (0.67), `rounding-1/3-at-2dp`
  (0.33), `rounding-5/8-at-3dp` (0.625), `rounding-zero-at-2dp`,
  `rounding-one-at-2dp`.
- `half-ties-round-down`, `half-ties-round-up` — the tie, in both directions.
- `agreement-2/3`, `agreement-1/3`, `agreement-5/8`, `agreement-one` — the two
  rules agree wherever there is no tie, so the tie direction is the *only* thing
  the contract decision can change.

These vectors are exported to `test/fixtures/agda-known-answers.json` for the
Julia conformance testset to reproduce.

## Open obligations

Everything not proved is in `proofs/residue/`:

- `exact-counts.residue` — `Checked.refused` injectivity (open, low impact); the
  `k = 63` ↔ Julia `Int64` correspondence (open by construction, closed by a
  Julia-side conformance assertion); the `:widen` branch (out of scope).
- `benjamini-hochberg.residue` — non-negativity of the scaled value (open, low
  impact, mechanical obstacle recorded); the step-down `min` envelope (not
  attempted; `ℚᵘ` has no `⊓` lemmas); FDR control (out of scope).
- `out-of-scope.residue` — probability theory, IEEE-754 semantics, the
  Agda-model-to-Julia correspondence, the `:ordinary` regression path, and exact
  statistical tests (issue #3).
- `toolchain.residue` — `R-TC-1`: the proofs compile only against an untagged
  agda-stdlib development SHA, not against any release; porting to `main`'s
  `v2.1` pin is required to land there, with the failing APIs listed in order.
  `R-TC-2` (closed): Agda's library-file location differs between install
  methods; `bootstrap.sh` now writes both locations and passes `--library-file`
  explicitly.

`q_i ≥ p_i` is **false** in general and is deliberately absent.

## Reproducing

```sh
proofs/bootstrap.sh            # install Agda + stdlib if absent, then check
just proofs                    # same, via the Justfile
just proofs-selftest           # prove the gate rejects broken proofs
```

The bootstrap pins Agda `2.7.0.1` and agda-stdlib at the SHA above. Changing
either is a change to the gate and must be reviewed — and note that moving the
stdlib pin to a *released* tag is not a version bump, it is a port (see `R-TC-1`).

`proofs/bootstrap.sh` fails non-zero if Agda cannot be installed. **An absent
prover is a failure, never a skip** — in CI and locally alike.
