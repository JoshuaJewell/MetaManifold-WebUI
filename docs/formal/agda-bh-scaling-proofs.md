<!--
SPDX-License-Identifier: AGPL-3.0-only
SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
-->

# Machine-checked proofs: Benjamini–Hochberg and size-factor scaling

`proofs/agda/` holds Agda proofs of the exact-arithmetic properties of three
functions in `src/analysis/differential.jl`: `bh_adjust`, `tss_factors` and
`rle_factors` (with `_median`). The functions were added by the
differential-abundance change (branch `feat/differential-abundance-nb`,
commits `6c3e8b2` and `23dc89c`). These proofs depend on that change being
merged first.

Every module is checked with `--safe --without-K`. There are no postulates,
no holes and no termination or positivity pragmas.
`proofs/tests/axiom-audit.sh` enforces this. It also fails if any module is
not imported by `MetaManifold/All.agda`, because a module that nothing imports
is never type-checked.

## Toolchain

The checker is **Agda 2.6.4.3** with **agda-stdlib 2.1**. The CI job
"Agda proofs" installs this pair from Debian trixie (`agda-bin=2.6.4.3-1+b2`,
`agda-stdlib=2.1-4`). It runs in a Debian image pinned by digest, because
Agda publishes no 2.6.4.3 binary and stdlib 2.1 needs Agda 2.6.4.x.

Earlier planning material and the archived proofs used Agda 2.7.0.1. These
proofs are written for, and checked with, 2.6.4.3. They have not been checked
with 2.7.0.1 or with stdlib 2.2 or later.

To check locally:

```sh
bash proofs/tests/axiom-audit.sh
cd proofs/agda && agda MetaManifold/All.agda   # needs standard-library-2.1 in $AGDA_DIR/libraries
```

## What is modelled

Values are exact rationals (`ℚ`), not `Float64`. A list stands in for a
vector.

- **`bhAdjust`** follows the shape of the Julia function:
  1. Sort the p-values, keeping each one's input index.
  2. Scale the value at rank `r` by `n / r`.
  3. Take the running minimum from the largest rank down, then clamp it at 1.
  4. Write each value back to its input position (`restore`).

  It differs from the Julia code in two ways:
  - **Running minimum.** Julia starts the running minimum at `1.0`; the model
    starts from the first scaled value and clamps afterwards.
    `stepUp-envelope` proves the two produce the same values.
  - **Tie order.** Julia sorts in descending order. Among equal p-values, its
    stable sort gives the earlier input the higher rank; the model's
    ascending `(p, index)` sort gives it the lower rank. Either way the
    sorted p-values are identical, and `bh-monotone` (applied in both
    directions) proves that equal p-values receive equal adjusted values. So
    the tie order cannot change the output. That last step is an argument
    about the Julia code built on a proved lemma. It is not a formal
    equivalence proof between Julia and Agda.
- **Input validation** (finite and in `[0, 1]`) is the hypothesis
  `All Valid ps`. The Julia `ArgumentError` paths are therefore assumptions
  of the theorems, not theorems.

## Theorem map

| Requirement | Agda theorem (module) | Plain statement | Julia function | D1 test (`test/unit/test_differential.jl`) |
|---|---|---|---|---|
| R-BH-1 | `bh-nonneg` (BenjaminiHochberg) | For valid p-values, every adjusted value is ≥ 0. | `bh_adjust` | lines 106–110, parity 117–124 |
| R-BH-2a | `bh-dominates-p` | Every adjusted value is ≥ its own p-value. | `bh_adjust` | line 107 (`0.001 → 0.003`), parity |
| R-BH-2b | `bh-at-most-one` | Every adjusted value is ≤ 1. | `bh_adjust` | line 107 (`0.5 → 0.5`), line 110 (`[1.0]`) |
| R-BH-2c | `bh-monotone` | For any two inputs, `p_i ≤ p_j` implies `adj_i ≤ adj_j`. This covers ties in both directions. | `bh_adjust` | line 109 (`[0.04, 0.01] → [0.04, 0.02]`) |
| R-BH-2d | `bh-envelope`, `stepUp-envelope` | In sorted order, the adjusted value at rank `k` is `min(1, min_{j ≥ k} p_(j) · n/j)`. | `bh_adjust` loop | line 106 (all → 0.04), parity |
| – | `bh-length` | The output has one entry per input. | `bh_adjust` | line 111 (empty → empty) |
| – | `monotone-in-family-size` | A larger family never gives a smaller scaled value. | `n / rank * p` | – |
| TSS | `total-positive` (Scaling) | A row of non-negative counts with at least one positive count has a positive total. So the `l > 0` refusal fires exactly on all-zero rows. | `tss_factors` | lines 137–139 (empty row refused) |
| TSS | `proportions-sum-to-1` | A row's proportions `x / Σx` sum to 1. | (relative abundance) | – |
| TSS | `proportions-bounded` | Every proportion lies in `[0, 1]`. | (relative abundance) | – |
| TSS/RLE | `centre-positive`, `centre-monotone` | Dividing factors by any positive `g` keeps them positive and keeps their order. | `lib ./ _geomean(lib)`, `raw ./ _geomean(raw)` | – |
| RLE | `ratio-positive` | A positive count over a positive normaliser is positive. | `view(sub, i, :) ./ gm` | – |
| RLE | `median-positive` | The median (taken as `_median` takes it) of a non-empty list of positive values is positive. | `_median`, `rle_factors` | lines 143–153 |

Taken together, the RLE rows show that once one taxon has reads in every
sample, every RLE factor is positive in exact arithmetic. This holds given a
positive geometric mean, which is assumed (see below).

`tss_factors` never forms proportions; it centres library sizes. The
proportion theorems are the relative-abundance reading of TSS, proved because
the requirement names them. They do not describe a value the Julia code
returns.

## Corrections to the requirement wording

- **R-BH-2a is true for the procedure.** An earlier archived note recorded
  "q_i ≥ p_i is FALSE". It is false only for a scaling factor below 1. In BH
  the rank `r` is at most `n`, so `n / r ≥ 1` (`rankFactor-≥1`), and the
  clamp at 1 keeps the bound because `p ≤ 1`. The theorem is proved as
  worded, for valid inputs.
- **R-BH-2c is stated over input positions, not over sorted ranks.** A
  sorted-rank statement would leave open how tied inputs are written back.
  The membership form `x ∈ zip ps (bhAdjust ps)` covers ties.
- **Empty input** is covered. All the `All`/membership statements are vacuous
  on `[]`, and `bh-length` gives the empty result.

## What is not proved, and why

- **Geometric means, `exp` and `log`.** A geometric mean of rationals is in
  general irrational, so `_geomean` cannot be computed over `ℚ`. Wherever
  Julia divides by a geometric mean, the theorems divide by an arbitrary
  positive rational. That the real geometric mean of positive numbers is
  positive is a fact about the reals, not a theorem here.
- **Scale invariance** (`tss_factors(7.5 .* x) ≈ tss_factors(x)`, line 136;
  the same for RLE, line 148). This is a property of geometric-mean
  centring and is not proved.
- **The RLE "doubled sample" ratio test** (lines 150–151) is not proved.
- **Median equivariance** (`median(c·x) = c·median(x)`, translation) is not
  proved.
- **Floating point.** Every statement is exact. `Float64` rounding, and a
  refusal caused by underflow to `0.0`, are outside the model.
- **Equivalence of the Agda and Julia code** is argued, not formalised. The
  model mirrors the loop, and the two differences above are each
  discharged by a proved lemma plus a short argument.

## Positive controls (run when the proofs were written)

1. **Audit.** Each of the following was planted in a scratch copy, and each
   made the audit exit 1 with a named finding: a `postulate` block, a
   mid-line `private postulate`, `?`, `{! !}`, `{-# TERMINATING #-}`, a
   module not imported by `All.agda`, and a removed `--without-K`. A comment
   that mentions "postulate" and "?" does not trigger it.
2. **Monotonicity needs the cumulative minimum.** Removing it makes Agda
   reject the proofs:
   - from `running`: rejected in `head-eq`, the envelope proof;
   - from the sorted stage (`Z`): rejected in `Z-bounded`;
   - from the stage the monotonicity proof alone uses: rejected in
     `Z-linked`, at `ℚₚ.p⊓q≤q s F`.

   A concrete counterexample is `p = [0.01, 0.011]`, `n = 2`. Without the
   minimum the scaled values are `0.02` and `0.011`, which are not monotone.
