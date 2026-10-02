-- SPDX-License-Identifier: AGPL-3.0-only
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- Library-size scaling (TSS) and median-of-ratios (RLE) over exact rationals:
-- the parts of `src/analysis/scaling.jl` that are statements about rational
-- arithmetic, proved; the parts that are not, named.
--
-- What is NOT modelled, and why:
--
--   * Geometric means, `exp` and `log`.  Julia centres every factor vector by
--     its geometric mean and turns factors into offsets with `log`; RLE divides
--     each count by the feature's geometric mean.  None of these is a rational
--     function, so none can be computed over ℚ.  Wherever Julia divides by a
--     geometric mean, the theorems below take the divisor as an ARBITRARY
--     positive rational `g`.  Every result therefore holds for the real
--     geometric mean as soon as it is positive (it is: a geometric mean of
--     positive numbers is positive), but that last step is a fact about the
--     reals, not a theorem here.
--   * Floating point.  Julia computes in Float64; these are exact statements.
--     `_require_positive` can fire in floating point (underflow to 0.0) where
--     the exact quantity is positive; it is not modelled.
--
-- Theorems (names as in docs/formal/agda-bh-scaling-proofs.md):
--
--   total-positive        a column of non-negative counts with at least one
--                         positive count has a positive total (so TSS refuses
--                         exactly the all-zero columns)
--   proportions-sum-to-1  the relative-abundance transform sums to 1
--   proportions-bounded   every proportion lies in [0, 1]
--   centre-positive       dividing by any positive g keeps every factor positive
--   centre-monotone       ... and keeps their order
--   ratio-positive        a positive count over a positive normaliser is positive
--   median-positive       the median of a non-empty list of positive ratios is
--                         positive (so RLE's `_require_positive` cannot fire in
--                         exact arithmetic once a usable feature exists)

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

-- The column total: Julia's `sum(counts[:, j])`.
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

-- A column of non-negative counts with at least one positive count has a
-- positive total.  So `_require_positive(totals, "tss", …)` refuses exactly
-- the columns whose counts are all zero.
total-positive : ∀ xs → All (0ℚ ≤_) xs → Any (0ℚ <_) xs → 0ℚ < sum xs
total-positive (x ∷ xs) (_   ∷ nn) (here 0<x)  =
  subst (_< x + sum xs) 0+0 (ℚₚ.+-mono-<-≤ 0<x (sum-nonNeg xs nn))
total-positive (x ∷ xs) (0≤x ∷ nn) (there any) =
  subst (_< x + sum xs) 0+0 (ℚₚ.+-mono-≤-< 0≤x (total-positive xs nn any))

-- The relative-abundance (TSS proportion) transform of one column.
proportions : (xs : List ℚ) → .{{Positive (sum xs)}} → List ℚ
proportions xs = map (_* inv (sum xs)) xs

-- The proportions of a column sum to one.
proportions-sum-to-1 : ∀ xs .{{_ : Positive (sum xs)}} →
  sum (proportions xs) ≡ 1ℚ
proportions-sum-to-1 xs =
  trans (sum-scale (inv (sum xs)) xs)
        (ℚₚ.*-inverseʳ (sum xs) {{ℚₚ.pos⇒nonZero (sum xs)}})

-- Every proportion of a column of non-negative counts lies in [0, 1].
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

-- The middle of a sorted list, as `Statistics.median` takes it: the element
-- at (one-based) position (n+1)/2 for odd n, the mean of the elements at
-- n/2 and n/2+1 for even n.  Both are at zero-based k = ⌊(n-1)/2⌋ (and k+1).
-- The defaults (the list's own head) are never read for a non-empty list.
middle : List ℚ → ℚ
middle []         = 0ℚ
middle (y ∷ ys)   =
  if isOdd n then at y k (y ∷ ys)
             else (at y k (y ∷ ys) + at y (suc k) (y ∷ ys)) * ½
  where
  n = suc (length ys)
  k = ⌊ length ys /2⌋   -- ⌊(n-1)/2⌋

-- Julia's `Statistics.median`: sort, then take the middle.
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

-- The median of a non-empty list of positive ratios is positive.  Hence RLE's
-- `_require_positive(ratios, …)` cannot fire in exact arithmetic once at
-- least one feature is positive in every sample (Julia refuses the case with
-- no such feature separately, before computing any median).
median-positive : ∀ xs → 0 ℕ.< length xs → All (0ℚ <_) xs → 0ℚ < median xs
median-positive xs 0<n ps =
  middle-positive (sort xs)
    (subst (0 ℕ.<_) (sym (↭ₚ.↭-length (sort-↭ xs))) 0<n)
    (↭ₚ.All-resp-↭ (↭-sym (sort-↭ xs)) ps)
