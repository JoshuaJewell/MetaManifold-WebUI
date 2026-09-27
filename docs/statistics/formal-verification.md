<!--
SPDX-License-Identifier: CC-BY-SA-4.0
-->

# Formal verification of the validated statistics layer

Status: **the Agda gate is green**. See [`proofs/PROOF-STATUS.md`](../../proofs/PROOF-STATUS.md)
for the full inventory and [`proofs/residue/`](../../proofs/residue/) for what is
explicitly *not* proved.

This page is for a reviewer who needs to know what "validated with a proof
assistant" does and does not mean for issue #1, without reading Agda.

## The one-line version

Seven modules, 77 top-level definitions, type-checked by Agda 2.7.0.1 with
agda-stdlib (pinned by SHA — see the caveat below) under `--safe --without-K`: no
postulates, no foreign code, no proof-irrelevance escape hatch, no universe-level
cheating. A separate audit script enforces those four facts on every build, and a
self-test script proves the audit and the type-checker can actually fail.

One caveat that changes what "verified" means here, stated up front rather than
in a footnote: the stdlib pin is `2ffa8b7d`, which is on the development line
towards agda-stdlib 3.0. **No `v3.0` tag exists** — the newest is `v2.4` — and
these proofs do not compile against `v2.4` or against the `v2.1` estate pin on
`main`. Reproducible, but not portable. `proofs/residue/toolchain.residue` lists
the exact APIs that have to change to port it.

## Why a proof assistant at all

Issue #1 asks for a numeric layer that never silently produces a plausible wrong
number. Tests sample the input space; they cannot cover it. Three of the claims
this layer makes are exactly the kind that survive every test anybody thinks to
write and fail in production:

1. **"An exact proportion of an empty sample is refused, not zero."** A test can
   check `c = 0, t = 0`. It cannot check that no future refactor turns the
   refusal into a default. In the model, the refusal is one arm of a sum type and
   `outcome-total` says there is no third arm — an edit that forgets a case does
   not type-check.

2. **"A count sum that overflows is refused, not wrapped."** `checkedAdd-is-exact`
   and `checkedAdd-refuses-exactly-when-it-must` are a two-sided statement: a
   returned value *is* the true sum, and a refusal happens *exactly* when the sum
   leaves the range. Proving one direction alone would be worthless, because each
   direction alone is satisfied by an implementation that is useless.

3. **"A permutation p-value cannot be zero."** `never-reports-zero` holds for
   *every* `b` and `B`. No rounding rule, display layer or downstream filter has
   to be trusted, because the estimator cannot reach zero. The negative control
   is proved alongside it: `naive-estimator-can-report-zero` shows the estimator
   this one replaces really does return exactly zero. A proof that a hazard was
   avoided is only evidence if the hazard is shown to be real.

## The design decision that made this tractable

`ℚᵘ` — Agda's unnormalised rationals — stores its denominator as
`suc denominator-1`. **A zero denominator is therefore unrepresentable by type.**

Consequences that fall out rather than being proved:

- `whole-is-one : ∀ n → mkℚᵘ (+[1+ n ]) n ≃ 1ℚᵘ` needs **no** `n ≢ 0` hypothesis.
- `value-has-nonzero-denominator` is a statement about a value that could not
  have been constructed otherwise.
- The refusal for division by zero is not a checked branch that could be
  reordered by a refactor; it is an impossibility in the representation.

This is the exact sense in which the refusals are total, and it is why Agda was
used rather than the Lean fallback the brief permitted.

## What a finding looked like

`relativeAbundance (suc c) 0` refuses as **`countExceedsTotal`**, not
`zeroTotal`. That is not a modelling choice: `exact_relative_abundance` checks
`count <= total` *before* `iszero(total)`, so a non-zero count against a zero
total reports the count error first. Both orders are defensible; only one is
implemented. `zero-total-with-count-is-refused-as-exceeding` records which, and a
test written from the docstring would have asserted the other.

## The bridge to Julia

A proof about a model proves nothing about an implementation unless something
ties them together. The tie is
[`test/fixtures/agda-known-answers.json`](../../test/fixtures/agda-known-answers.json):
values the proof assistant has *checked against a definition*, each annotated
with the lemma that fixes it. The Julia conformance testset must reproduce them
exactly. If the two disagree, the Julia layer is wrong.

The rounding vectors are the clearest case. "Correctly rounded to *s* decimals"
is specified as a predicate over integers alone — no division, so no rounding
inside the specification of rounding — and each printed decimal is then a
type-checked instance of it. `0.67` for `2/3` at two decimals is not a number
somebody typed into two places.

The same file records the tie-handling contract: at an exact half,
`IsRounding` goes down (`1/2` at 0 dp is `0`) and `IsRoundHalfUp` goes up
(`1/2` at 0 dp is `1`). The two agree everywhere else, which is itself proved, so
the tie direction is the *only* thing the contract decision can change.

## What this does not prove

Stated here because the alternative is a reader inferring it from the presence of
a `proofs/` directory:

- **No probability theory.** Nothing about coverage, bias, type I error or FDR
  control. Those are simulation-study questions with pre-set tolerances, and they
  live in the Julia validation layer.
- **No IEEE-754.** Every theorem is over exact rationals. No claim of the form
  "Float64 gives the same answer" is made anywhere.
- **No proof that `numeric_policy.jl` implements the model.** The bridge is the
  fixture file and the conformance testset — a test, not a proof.
- **Nothing about the `:ordinary` path.** Unchanged existing behaviour is a
  regression property over the existing suite.
- **`q_i ≥ p_i` is false** in general and is deliberately absent.

Each of these has an entry in `proofs/residue/` with an id, a precise statement,
and what would close it.

## Running it

```sh
just proofs          # bootstrap if needed, audit, type-check everything
just proofs-audit    # escape-hatch audit only
just proofs-selftest # prove the gate rejects broken proofs (nine breakages)
just proofs-clean    # remove the vendored toolchain and interface cache
```

`proofs/bootstrap.sh` pins Agda `2.7.0.1` and agda-stdlib at `2ffa8b7d`, and
exits non-zero if it cannot install them. **An absent prover is a failure, never
a skip**, in CI and locally alike. CI runs the audit, the type-check, the
self-test, and a consistency check that `PROOF-STATUS.md` matches the tree, and
uploads the type-check transcript as an artifact whether it passed or not.

Verified from an empty `proofs/.vendor`, so the numbers above are the bootstrap
script installing its own toolchain, not a pre-existing one being reused.
