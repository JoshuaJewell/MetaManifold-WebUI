<!-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk> -->
<!-- SPDX-License-Identifier: CC-BY-SA-4.0 -->

# Descriptive statistics contracts

This opt-in layer summarises a features-by-samples table. It does not perform
hypothesis tests, confidence intervals, differential abundance, normalisation,
or correction for sampling and measurement error. Existing R-backed analyses
continue to own their statistical methods.

For non-negative integer counts `c[i,j]`, the sample total is
`T[j] = sum_i c[i,j]`, evaluated with `BigInt`. If `T[j] > 0`, each proportion
is the reduced `Rational{BigInt}` value `c[i,j] / T[j]`. Their sum is exactly
one and multiplying any proportion by its total reconstructs its count.
For `T[j] = 0`, proportions are undefined (`null`), including when every count
is zero. Observed richness counts strictly positive entries. Groups pool
counts before division; they do not average sample proportions.

High-precision mode computes Shannon entropy `-sum(p * log(p))` with natural
logarithms and Gini–Simpson diversity `1 - sum(p^2)`. The latter uses an exact
rational intermediate before its final MPFR conversion. Both returned
diversity values are labelled approximations. Arithmetic precision is not a
confidence level. Zero-total diversity values are undefined.

The implementation scopes both MPFR precision and RoundNearest, restoring
the caller's settings even on an exception. A lock serialises this layer's
precision changes. It cannot stop unrelated code outside the lock from changing
process MPFR defaults. The compiler, GMP, MPFR, Julia's rational implementation,
and hardware remain part of the trusted execution base.

## Boundaries

* Direct HTTP input must use integer JSON tokens through `2^53 - 1`; larger
  counts must be whole-number strings. Fractional/exponent tokens, including
  `1.0` and `1e0`, are refused before their syntax can be erased by JSON parsing.
* The Julia exact-summary API can explicitly accept integral floating sources
  for compatibility, but labels the source approximate. HTTP and biological
  validation use `allow_approximate=false`.
* DuckDB input must use integer columns. Addition and grouping use Julia
  `BigInt`, not SQL `SUM` or a Float64 matrix. NULL and negative counts fail.
  Numeric filters on count columns and sample read-total filters fail.
* The run route applies explicitly requested supported row filters. It does
  not inherit category exclusions or diversity normalisation. Feature-level
  counts are the default; an explicit taxonomy aggregation describes those
  labels and does not establish their biological correctness.
* Counts outside the browser-safe range and all rational values are strings.
  The React view never converts them through JavaScript Number. Decimal
  displays are explicitly rounded; exact rational values remain available.
* Derived records use canonical JSON, content hashes and exclusive atomic
  publication. Existing conflicting or corrupt records are refused. Source
  result databases are opened read-only.
* Cancellation is sticky. Time limits are cooperative checks between work
  units, not hard interrupts inside GMP, MPFR, DuckDB or a filesystem call.

The executable limits are in `NumericPolicy` and the method catalogue tests
require matching advertised limits: 100,000 cells, 16 MiB request bodies,
4,096-bit input counts and MPFR precision, 65,536-bit rational intermediates,
100 display decimal places, and 64 MiB persisted records. A limit failure
returns an explicit error; it never substitutes a floating method.

## Mathematical evidence and its limits

The exact operations have the following algebraic obligations, conditional on
the validated nonnegative integer inputs and a positive total `T`:

* Proportions: `sum_i (c[i]/T) = (sum_i c[i])/T = 1`, and
  `T * (c[i]/T) = c[i]`. Reducing each fraction does not change its value.
* Pooling a group `G`: `C[i] = sum_(j in G) c[i,j]` and
  `sum_i C[i] = sum_(j in G) T[j]`, by exchanging finite integer sums.
  Dividing these pooled counts gives the pooled proportions; averaging the
  original proportions is generally a different operation.
* The exact Gini–Simpson intermediate is
  `(T^2 - sum_i c[i]^2)/T^2 = 2*sum_(i<k) c[i]*c[k]/T^2`.
  Its numerator is nonnegative and at most `T^2`, so this rational lies in
  `[0,1]`. Its final MPFR conversion is still labelled approximate.
* For display at `k` decimal places, set `n = abs(numerator)*10^k` and let
  `d` be the positive denominator. Divide exactly: `n = q*d + r`, with
  `0 <= r < d`. Choosing `q` when `2*r < d` and `q+1` otherwise gives a
  displayed value with absolute error at most `10^(-k)/2`.
  Restoring the original sign makes an exact tie round away
  from zero. This is rounding of the represented value, not an error bound
  relative to an unknown true value behind a floating-point measurement.

These are proofs of the stated arithmetic identities and rounding rule. They
do not remove the need to establish that parsing, grouping, computation and
display actually implement the corresponding operations. For `T = 0`, the
division obligations do not apply; the API must report undefined values.

The related [ExactCounts.jl proof model](https://github.com/hyperpolymath/ExactCounts.jl/blob/08c5c641775bf96508d1ee845e580c88c6d8e29a/proofs/PROOF-STATUS.md)
contains Agda proofs of checked addition,
proportions summing to one, aggregation over a common denominator, and
round-half-away-from-zero. Its proof-status document distinguishes the model
theorems from Julia tests. This application does not add an unregistered Julia
dependency or vendor that package. Source revisions and observed checker
results belong in the accompanying validation receipt.

Those model theorems are not an end-to-end verification of this application.
Julia property tests, independent Python `Fraction` and `Decimal` checks,
browser boundary tests and DuckDB integration tests test the implementation
connection. They do not prove an absence of all software defects, statistical
calibration, biological representativeness, or permission to redistribute data.
