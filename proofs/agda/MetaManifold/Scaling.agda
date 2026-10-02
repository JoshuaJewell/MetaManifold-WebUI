-- SPDX-License-Identifier: AGPL-3.0-only
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- Library-size scaling (TSS) and median-of-ratios (RLE) over exact rationals:
-- the parts of `tss_factors`, `rle_factors` and `_median` in
-- `src/analysis/differential.jl` (branch feat/differential-abundance-nb,
-- commit 6c3e8b2) that are statements about rational arithmetic, proved; the
-- parts that are not, named.
--
-- What is NOT modelled, and why:
--
--   * Geometric means, `exp` and `log`.  Julia centres each factor vector by
--     its geometric mean (`_geomean`), RLE divides each count by its taxon's
--     geometric mean, and the model uses `log` of the factors as offsets.
--     None of these is a rational function, and a geometric mean of
--     rationals is in general irrational, so none can be computed over ℚ.
--     Wherever Julia divides by a geometric mean, the theorems below divide
--     by an ARBITRARY positive rational `g`.  The same algebra holds over the
--     reals for any positive divisor, but that transfer is an argument, not a
--     theorem here.  In particular the scale invariance D1 tests
--     (`tss_factors(7.5 .* fixture) ≈ f`) is a property of the geometric
--     mean and is not proved.
--   * Floating point.  Julia computes in Float64; these are exact statements.
--     A refusal can fire in floating point (underflow to 0.0) where the exact
--     quantity is positive; that is not modelled.
--
-- D1's `tss_factors` computes library sizes (row sums of the samples x taxa
-- matrix) and centres them; it never forms proportions.  The proportion
-- theorems below are the relative-abundance reading of TSS asked for in the
-- D2 brief, and hold for any row with a positive total.
--
-- Theorems (names as in proofs/agda/README.md):
--
--   total-positive        a row of non-negative counts with at least one
--                         positive count has a positive total (so TSS's
--                         `l > 0` refusal fires exactly on all-zero rows)
--   proportions-sum-to-1  a row's proportions sum to 1
--   proportions-bounded   every proportion lies in [0, 1]
--   centre-positive       dividing by any positive g keeps every factor positive
--   centre-monotone       ... and keeps their order
--   ratio-positive        a positive count over a positive normaliser is positive
--   median-positive       the median (as `_median` takes it) of a non-empty list
--                         of positive ratios is positive, so once one taxon has
--                         reads in every sample, every RLE factor is positive in
--                         exact arithmetic

{-# OPTIONS --safe --without-K #-}

module MetaManifold.Scaling where

open import Data.Product.Base using (_×_; _,_)
open import Data.Bool.Base using (Bool; true; false; not; if_then_else_)
open import Data.Nat.Base as ℕ using (ℕ; zero; suc; ⌊_/2⌋)
open import Data.List.Base using (List; []; _∷_; map; foldr; length)
open import Data.List.Relation.Unary.All as All using (All; []; _∷_)
open import Data.List.Relation.Unary.Any using (Any; here; there)
import Data.List.Relation.Binary.Permutation.Propositional.Properties as ↭ₚ
open import Data.List.Relation.Binary.Permutation.Propositional using (↭-sym)
open import Data.Rational.Base as ℚ
  using (ℚ; 0ℚ; 1ℚ; ½; _+_; _*_; 1/_; _≤_; _<_; Positive)
import Data.Rational.Properties as ℚₚ
open import Relation.Binary.PropositionalEquality
  using (_≡_; refl; sym; trans; cong; subst)

open import Data.List.Sort.MergeSort ℚₚ.≤-decTotalOrder using (sort; sort-↭)

------------------------------------------------------------------------
-- Shared arithmetic

-- A sample's library size: Julia's `vec(sum(counts; dims=2))`, one row.
sum : List ℚ → ℚ
sum = foldr _+_ 0ℚ

-- Division by a positive rational, as multiplication by its inverse.
inv : (g : ℚ) → .{{Positive g}} → ℚ
inv g = (1/ g) {{ℚₚ.pos⇒nonZero g}}

private
  inv-positive : ∀ g .{{_ : Positive g}} → Positive (inv g)
  inv-positive g = ℚₚ.1/pos⇒pos g

  inv-nonNeg : ∀ g .{{_ : Positive g}} → 0ℚ ≤ inv g
  inv-nonNeg g = ℚₚ.<⇒≤ (ℚₚ.positive⁻¹ (inv g) {{inv-positive g}})

  -- x ↦ x · r is monotone for a non-negative r.
  scale-≤ : ∀ r → 0ℚ ≤ r → ∀ {x y} → x ≤ y → x * r ≤ y * r
  scale-≤ r 0≤r = ℚₚ.*-monoʳ-≤-nonNeg r {{ℚ.nonNegative 0≤r}}

  -- x ↦ x · r is strictly monotone for a positive r.
  scale-< : ∀ r → .{{_ : Positive r}} → ∀ {x y} → x < y → x * r < y * r
  scale-< r = ℚₚ.*-monoˡ-<-pos r

  0+0 : 0ℚ + 0ℚ ≡ 0ℚ
  0+0 = ℚₚ.+-identityʳ 0ℚ

------------------------------------------------------------------------
-- TSS: totals and proportions

private
  sum-nonNeg : ∀ xs → All (0ℚ ≤_) xs → 0ℚ ≤ sum xs
  sum-nonNeg []       []         = ℚₚ.≤-refl
  sum-nonNeg (x ∷ xs) (0≤x ∷ nn) =
    subst (_≤ x + sum xs) 0+0 (ℚₚ.+-mono-≤ 0≤x (sum-nonNeg xs nn))

  -- Every non-negative summand is at most the total.
  below-sum : ∀ xs → All (0ℚ ≤_) xs → All (_≤ sum xs) xs
  below-sum []       []         = []
  below-sum (x ∷ xs) (0≤x ∷ nn) =
    ℚₚ.≤-trans (ℚₚ.≤-reflexive (sym (ℚₚ.+-identityʳ x)))
               (ℚₚ.+-mono-≤ (ℚₚ.≤-refl {x}) (sum-nonNeg xs nn))
    ∷ All.map (λ y≤s → ℚₚ.≤-trans y≤s
                 (ℚₚ.≤-trans (ℚₚ.≤-reflexive (sym (ℚₚ.+-identityˡ (sum xs))))
                             (ℚₚ.+-mono-≤ 0≤x (ℚₚ.≤-refl {sum xs}))))
              (below-sum xs nn)

  sum-scale : ∀ r xs → sum (map (_* r) xs) ≡ sum xs * r
  sum-scale r []       = sym (ℚₚ.*-zeroˡ r)
  sum-scale r (x ∷ xs) =
    trans (cong (x * r +_) (sum-scale r xs)) (sym (ℚₚ.*-distribʳ-+ r x (sum xs)))

-- A row of non-negative counts with at least one positive count has a
-- positive total.  So `tss_factors`'s `l > 0` refusal fires exactly on the
-- rows whose counts are all zero (`validate_counts` guarantees non-negative).
total-positive : ∀ xs → All (0ℚ ≤_) xs → Any (0ℚ <_) xs → 0ℚ < sum xs
total-positive (x ∷ xs) (_   ∷ nn) (here 0<x)  =
  subst (_< x + sum xs) 0+0 (ℚₚ.+-mono-<-≤ 0<x (sum-nonNeg xs nn))
total-positive (x ∷ xs) (0≤x ∷ nn) (there any) =
  subst (_< x + sum xs) 0+0 (ℚₚ.+-mono-≤-< 0≤x (total-positive xs nn any))

-- The relative-abundance (TSS proportion) transform of one row.
proportions : (xs : List ℚ) → .{{Positive (sum xs)}} → List ℚ
proportions xs = map (_* inv (sum xs)) xs

-- The proportions of a row sum to one.
proportions-sum-to-1 : ∀ xs .{{_ : Positive (sum xs)}} →
  sum (proportions xs) ≡ 1ℚ
proportions-sum-to-1 xs =
  trans (sum-scale (inv (sum xs)) xs)
        (ℚₚ.*-inverseʳ (sum xs) {{ℚₚ.pos⇒nonZero (sum xs)}})

-- Every proportion of a row of non-negative counts lies in [0, 1].
proportions-bounded : ∀ xs .{{_ : Positive (sum xs)}} → All (0ℚ ≤_) xs →
  All (λ p → 0ℚ ≤ p × p ≤ 1ℚ) (proportions xs)
proportions-bounded xs nn = go xs nn (below-sum xs nn)
  where
  T = sum xs
  r = inv T
  0≤r = inv-nonNeg T
  go : ∀ ys → All (0ℚ ≤_) ys → All (_≤ T) ys →
       All (λ p → 0ℚ ≤ p × p ≤ 1ℚ) (map (_* r) ys)
  go []       []         []         = []
  go (y ∷ ys) (0≤y ∷ nn′) (y≤T ∷ bs) =
    ( ℚₚ.≤-trans (ℚₚ.≤-reflexive (sym (ℚₚ.*-zeroˡ r))) (scale-≤ r 0≤r 0≤y)
    , ℚₚ.≤-trans (scale-≤ r 0≤r y≤T)
                 (ℚₚ.≤-reflexive (ℚₚ.*-inverseʳ T {{ℚₚ.pos⇒nonZero T}})) )
    ∷ go ys nn′ bs

------------------------------------------------------------------------
-- Centring: Julia's `raw ./ g` for a positive centring constant g

-- Divide every raw quantity by the same positive constant.
centre : (g : ℚ) → .{{Positive g}} → List ℚ → List ℚ
centre g = map (_* inv g)

-- Centring keeps every factor positive.
centre-positive : ∀ g .{{_ : Positive g}} xs →
  All (0ℚ <_) xs → All (0ℚ <_) (centre g xs)
centre-positive g []       []         = []
centre-positive g (x ∷ xs) (0<x ∷ ps) =
  subst (_< x * inv g) (ℚₚ.*-zeroˡ (inv g))
        (scale-< (inv g) {{inv-positive g}} 0<x)
  ∷ centre-positive g xs ps

-- Centring keeps the order of any two factors.
centre-monotone : ∀ g .{{_ : Positive g}} {x y} → x ≤ y → x * inv g ≤ y * inv g
centre-monotone g = scale-≤ (inv g) (inv-nonNeg g)

------------------------------------------------------------------------
-- RLE: ratios and their median

-- One RLE ratio: a count over its feature's normaliser (the geometric mean in
-- Julia; any positive rational here).
ratio-positive : ∀ c g .{{_ : Positive g}} → 0ℚ < c → 0ℚ < c * inv g
ratio-positive c g 0<c =
  subst (_< c * inv g) (ℚₚ.*-zeroˡ (inv g))
        (scale-< (inv g) {{inv-positive g}} 0<c)

-- Parity, by recursion (so that proofs can follow the same case split).
isOdd : ℕ → Bool
isOdd zero    = false
isOdd (suc n) = not (isOdd n)

-- The element at zero-based position k, or the default d past the end.
at : ℚ → ℕ → List ℚ → ℚ
at d k       []       = d
at d zero    (y ∷ ys) = y
at d (suc k) (y ∷ ys) = at d k ys

-- The middle of a sorted list, as D1's `_median` takes it: with m = n ÷ 2,
-- `s[m + 1]` for odd n and `(s[m] + s[m + 1]) / 2` for even n (one-based).
-- Both are at zero-based k = ⌊(n-1)/2⌋ (and k+1).  `at` falls back to the
-- list's own head past the end; for these k it is not reached (not proved
-- here; positivity does not depend on it).
middle : List ℚ → ℚ
middle []         = 0ℚ
middle (y ∷ ys)   =
  if isOdd n then at y k (y ∷ ys)
             else (at y k (y ∷ ys) + at y (suc k) (y ∷ ys)) * ½
  where
  n = suc (length ys)
  k = ⌊ length ys /2⌋   -- ⌊(n-1)/2⌋

-- D1's `_median`: sort, then take the middle.
median : List ℚ → ℚ
median xs = middle (sort xs)

private
  at-All : ∀ {P : ℚ → Set} d k ys → P d → All P ys → P (at d k ys)
  at-All d k       []       pd _          = pd
  at-All d zero    (y ∷ ys) pd (py ∷ _)   = py
  at-All d (suc k) (y ∷ ys) pd (_  ∷ pys) = at-All d k ys pd pys

  mean-positive : ∀ {a b} → 0ℚ < a → 0ℚ < b → 0ℚ < (a + b) * ½
  mean-positive {a} {b} 0<a 0<b =
    subst (_< (a + b) * ½) (ℚₚ.*-zeroˡ ½)
      (scale-< ½ (subst (_< a + b) 0+0 (ℚₚ.+-mono-< 0<a 0<b)))

  middle-positive : ∀ ys → 0 ℕ.< length ys → All (0ℚ <_) ys → 0ℚ < middle ys
  middle-positive (y ∷ ys) _ ps@(py ∷ _) with isOdd (suc (length ys))
  ... | true  = at-All y k (y ∷ ys) py ps
    where k = ⌊ length ys /2⌋
  ... | false = mean-positive (at-All y k (y ∷ ys) py ps)
                              (at-All y (suc k) (y ∷ ys) py ps)
    where k = ⌊ length ys /2⌋

-- The median of a non-empty list of positive ratios is positive.  Hence every RLE
-- factor is positive in exact arithmetic once one taxon has reads in every
-- sample (`rle_factors` refuses the case with no such taxon before computing
-- any median).
median-positive : ∀ xs → 0 ℕ.< length xs → All (0ℚ <_) xs → 0ℚ < median xs
median-positive xs 0<n ps =
  middle-positive (sort xs)
    (subst (0 ℕ.<_) (sym (↭ₚ.↭-length (sort-↭ xs))) 0<n)
    (↭ₚ.All-resp-↭ (↭-sym (sort-↭ xs)) ps)
