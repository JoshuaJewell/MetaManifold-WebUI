<!--
SPDX-License-Identifier: AGPL-3.0-only
SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
-->

# Machine-checked proofs: Benjamini–Hochberg, size-factor scaling and PERMANOVA permutations

`proofs/agda/` holds Agda proofs of the exact-arithmetic properties of three
functions in `src/analysis/differential.jl`: `bh_adjust`, `tss_factors` and
`rle_factors` (with `_median`).

`MetaManifold/Permanova.agda` covers the two things `run_permanova` in
`src/analysis/analysis.jl` decides for itself rather than taking from vegan:
the reported `min_p_value`, and which terms have their p-value withheld
(`_varies_within_any_block`).

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

The proofs have not been checked with Agda 2.7 or with stdlib 2.2 or later.

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
- **`variesWithinAnyBlock`** follows the shape of
  `_varies_within_any_block`. A list of `(value, block)` pairs stands in for
  `zip(values, blocks)` and an association list for the `first_seen`
  dictionary. As with `get!`, a block's value is stored only the first time
  the block is seen, and the loop returns `true` at the first sample whose
  value differs from its block's stored one. Equality of values and of blocks
  is any decidable equality, standing in for `isequal`.
- **Within-block permutations** are permutations of the samples
  (`Permutation′ n` from the standard library) that keep every sample in its
  block. A term's labels are a function from samples to values. The
  p-value is vegan's `(k + 1)/(B + 1)` for `B` permutations, `k` of which
  reached the observed statistic.

## Theorem map

### Benjamini–Hochberg and scaling (`BenjaminiHochberg`, `Scaling`)

Tests are testsets in `test/unit/test_differential.jl`.

| Agda theorem | Plain statement | Julia function | Test |
|---|---|---|---|
| `bh-nonneg` | For valid p-values, every adjusted value is ≥ 0. | `bh_adjust` | "BH adjustment: known answers", "parity with R p.adjust" |
| `bh-dominates-p` | Every adjusted value is ≥ its own p-value. | `bh_adjust` | "BH adjustment: known answers" (`0.001 → 0.003`), parity |
| `bh-at-most-one` | Every adjusted value is ≤ 1. | `bh_adjust` | "BH adjustment: known answers" (`0.5 → 0.5`, `[1.0]`) |
| `bh-monotone` | For any two inputs, `p_i ≤ p_j` implies `adj_i ≤ adj_j`. This covers ties in both directions. | `bh_adjust` | "BH adjustment: known answers" (`[0.04, 0.01] → [0.04, 0.02]`) |
| `bh-envelope`, `stepUp-envelope` | In sorted order, the adjusted value at rank `k` is `min(1, min_{j ≥ k} p_(j) · n/j)`. | `bh_adjust` loop | "BH adjustment: known answers" (all → 0.04), parity |
| `bh-length` | The output has one entry per input; empty input gives empty output. | `bh_adjust` | "BH adjustment: known answers" (empty → empty) |
| `monotone-in-family-size` | A larger family never gives a smaller scaled value. | `n / rank * p` | – |
| `total-positive` | A row of non-negative counts with at least one positive count has a positive total. So the `l > 0` refusal fires exactly on all-zero rows. | `tss_factors` | "TSS size factors" (empty row refused) |
| `proportions-sum-to-1` | A row's proportions `x / Σx` sum to 1. | (relative abundance) | – |
| `proportions-bounded` | Every proportion lies in `[0, 1]`. | (relative abundance) | – |
| `centre-positive`, `centre-monotone` | Dividing factors by any positive `g` keeps them positive and keeps their order. | `lib ./ _geomean(lib)`, `raw ./ _geomean(raw)` | – |
| `ratio-positive` | A positive count over a positive normaliser is positive. | `view(sub, i, :) ./ gm` | – |
| `median-positive` | The median (taken as `_median` takes it) of a non-empty list of positive values is positive. | `_median`, `rle_factors` | "RLE size factors" |

Taken together, the RLE rows show that once one taxon has reads in every
sample, every RLE factor is positive in exact arithmetic. This holds given a
positive geometric mean, which is assumed (see below).

`tss_factors` never forms proportions; it centres library sizes. The
proportion theorems are the relative-abundance reading of TSS and do not
describe a value the Julia code returns.

### PERMANOVA (`Permanova`)

Tests are testsets in `test/unit/test_permanova.jl`.

| Agda theorem | Plain statement | Julia | Test |
|---|---|---|---|
| `minP-attained`, `p-floor` | `1/(B + 1)` is the smallest p-value `B` permutations can give, and it is given when no permutation reaches the observed statistic. | `min_p_value = 1 / (Int(n_perm) + 1)` | "each term is tested on its own", "within-individual blocks: permutation count, floor and refusal" |
| `p-at-most-one` | With `k ≤ B`, no p-value is above 1. | (vegan) | – |
| `below-floor-unreachable` | A threshold below `1/(B + 1)` is below every p-value, so no term can meet it. | `min_p_value` | – |
| `minP-anti` | More permutations never raise the floor. | – | – |
| `minP-63`, `no-p-at-1%-with-63` | With 63 permutations the floor is `1/64`, and no p-value is at or below 0.01. | – | "within-individual blocks: permutation count, floor and refusal" |
| `never-exchanged` | A term that is constant within every block is left in place by every within-block permutation, so its labels are never exchanged. | the reason for withholding | – |
| `exchanged` | A term with two values in one block is moved by a within-block permutation: the swap of those two samples. | – | – |
| `untestable⇔fixed` | Constant within every block ⇔ fixed by every within-block permutation. | – | – |
| `varies-false⇔Constant` | The model of `_varies_within_any_block` returns `false` exactly when no block holds two different values. | `_varies_within_any_block` | "_varies_within_any_block" |
| `withheld⇔never-exchanged` | So `untestable` (the function returned `false`) holds exactly when no within-block permutation exchanges the term's labels. | `untestable`, `withheld` | "within-individual blocks: permutation count, floor and refusal" |
| `tested⇔exchanged` | The block-constancy filter admits a term exactly when some within-block permutation changes its labels; this does not guarantee a reported p-value. | `untestable`, `withheld` | "within-individual blocks: permutation count, floor and refusal" |
| `constant-in-blocks`, `varies-in-a-block` | The two direct cases in the test, by evaluation. | `_varies_within_any_block` | "_varies_within_any_block" |
| `fixture-group-withheld`, `fixture-run-tested` | On the test's twelve samples (six individuals, two runs each, three individuals per group), `group` is never exchanged and `run` is. | `run_permanova` | "within-individual blocks: permutation count, floor and refusal" |

## What is not proved, and why

- **Statistical validity and successful reporting.** Label exchangeability here
  characterises only the block-constancy filter. It does not prove model
  identifiability, assumptions for inference, or that vegan returns a p-value.
  Julia also withholds missing p-values independently of this filter.

- **Geometric means, `exp` and `log`.** A geometric mean of rationals is in
  general irrational, so `_geomean` cannot be computed over `ℚ`. Wherever
  Julia divides by a geometric mean, the theorems divide by an arbitrary
  positive rational. That the real geometric mean of positive numbers is
  positive is a fact about the reals, not a theorem here.
- **Scale invariance** (`tss_factors(7.5 .* x) ≈ tss_factors(x)`, and the
  same for RLE). This is a property of geometric-mean centring and is not
  proved.
- **The RLE "doubled sample" ratio test** is not proved.
- **Median equivariance** (`median(c·x) = c·median(x)`, translation) is not
  proved.
- **Floating point.** Every statement is exact. `Float64` rounding, and a
  refusal caused by underflow to `0.0`, are outside the model.
- **PERMANOVA: the pseudo-F statistic, and how many within-block
  permutations exist.** The p-value theorems take `B` and `k` as given. `B`
  is what `permute` reports, and `k` depends on the F distribution, neither of
  which is modelled. Six blocks of two samples admit `2⁶ = 64`
  arrangements, the identity and 63 others, which is the 63 in the test.
  That count (the product of the blocks' factorials) is not proved.
- **PERMANOVA: the unblocked case.** Julia asks `_varies_within_any_block`
  only when some block label repeats (`blocked`). Otherwise the permutations
  are unrestricted and nothing is withheld; that is Julia's `if`, not a
  theorem.
- **Equivalence of the Agda and Julia code** is argued, not formalised. The
  model mirrors the loop, and the two differences above are each
  discharged by a proved lemma plus a short argument.
