-- SPDX-License-Identifier: AGPL-3.0-only
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- Shared arithmetic for the machine-checked statistics core.
--
-- Two kinds of fact live here:
--
--   * `numerator-monotone-≤`, over the unnormalised rationals `ℚᵘ`, ported from
--     the archived proof tree because the Benjamini–Hochberg scaling lemmas are
--     stated over `ℚᵘ` and rest on it;
--   * the three facts about the Benjamini–Hochberg rank factor `m/(j+1)` that
--     every BH theorem needs: it is non-negative, it is at least one whenever
--     the rank is inside the family, and it shrinks as the rank grows.
--
-- The rank-factor facts are stated over the normalised rationals `ℚ` (the type
-- the BH procedure itself is defined over) and proved by moving to `ℚᵘ`, where
-- `m/(j+1)` is literally `mkℚᵘ (+ m) j` and a comparison is a single
-- cross-multiplication of integers.
--
-- Nothing in `proofs/agda` assumes anything beyond the Agda standard library.
-- `proofs/tests/axiom-audit.sh` checks that there is no postulate, hole or
-- unsafe pragma anywhere in the tree.

{-# OPTIONS --safe --without-K #-}

module MetaManifold.Prelude where

open import Data.Integer.Base as ℤ using (ℤ; +_; +[1+_])
import Data.Integer.Properties as ℤₚ
open import Data.Nat.Base as ℕ using (ℕ; zero; suc; z≤n; s≤s)
import Data.Nat.Properties as ℕₚ
open import Data.Rational.Base as ℚ
  using (ℚ; _/_; toℚᵘ; fromℚᵘ; 0ℚ; 1ℚ)
import Data.Rational.Properties as ℚₚ
open import Data.Rational.Unnormalised.Base as ℚᵘ
  using (ℚᵘ; mkℚᵘ; *≤*)
import Data.Rational.Unnormalised.Properties as ℚᵘₚ
open import Relation.Binary.PropositionalEquality using (_≡_; refl; sym; cong)

------------------------------------------------------------------------
-- Ported from the archive (over ℚᵘ)

-- Over a common denominator, a larger numerator is a larger fraction.
numerator-monotone-≤ :
  ∀ (n m : ℤ) (d : ℕ) → n ℤ.≤ m → mkℚᵘ n d ℚᵘ.≤ mkℚᵘ m d
numerator-monotone-≤ n m d n≤m =
  *≤* (ℤₚ.*-monoʳ-≤-nonNeg +[1+ d ] n≤m)

------------------------------------------------------------------------
-- The Benjamini–Hochberg rank factor m/(j+1)

-- The factor that multiplies the p-value at zero-based rank `j` in a family of
-- `m` tests.  Julia writes `n / i` with a one-based `i`; this is the same
-- number with `i = j + 1`.
rankFactor : ℕ → ℕ → ℚ
rankFactor m j = + m / suc j


-- Moving an inequality from ℚᵘ to ℚ through normalisation.
fromℚᵘ-mono-≤ : ∀ {p q} → p ℚᵘ.≤ q → fromℚᵘ p ℚ.≤ fromℚᵘ q
fromℚᵘ-mono-≤ {p} {q} p≤q = ℚₚ.toℚᵘ-cancel-≤
  (ℚᵘₚ.≤-respˡ-≃ (ℚᵘₚ.≃-sym (ℚₚ.toℚᵘ-fromℚᵘ p))
    (ℚᵘₚ.≤-respʳ-≃ (ℚᵘₚ.≃-sym (ℚₚ.toℚᵘ-fromℚᵘ q)) p≤q))

-- A rank factor is never negative.
rankFactor-nonNeg : ∀ m j → 0ℚ ℚ.≤ rankFactor m j
rankFactor-nonNeg m j = fromℚᵘ-mono-≤ {mkℚᵘ (+ 0) 0} {mkℚᵘ (+ m) j}
  (*≤* (ℤₚ.≤-trans (ℤₚ.≤-reflexive (ℤₚ.*-zeroˡ +[1+ j ]))
                   (ℤₚ.≤-trans (ℤ.+≤+ z≤n)
                     (ℤₚ.≤-reflexive (sym (ℤₚ.*-identityʳ (+ m)))))))

-- Inside the family (rank j+1 ≤ m) the rank factor is at least one, so
-- scaling never lowers a non-negative p-value.
rankFactor-≥1 : ∀ m j → suc j ℕ.≤ m → 1ℚ ℚ.≤ rankFactor m j
rankFactor-≥1 m j j<m = fromℚᵘ-mono-≤ {mkℚᵘ (+ 1) 0} {mkℚᵘ (+ m) j}
  (*≤* (ℤₚ.≤-trans (ℤₚ.≤-reflexive (ℤₚ.*-identityˡ +[1+ j ]))
         (ℤₚ.≤-trans (ℤ.+≤+ j<m)
                     (ℤₚ.≤-reflexive (sym (ℤₚ.*-identityʳ (+ m)))))))

-- A later rank has a smaller (or equal) factor: m/(j+2) ≤ m/(j+1).
rankFactor-anti : ∀ m j → rankFactor m (suc j) ℚ.≤ rankFactor m j
rankFactor-anti m j = fromℚᵘ-mono-≤ {mkℚᵘ (+ m) (suc j)} {mkℚᵘ (+ m) j}
  (*≤* (ℤₚ.*-monoˡ-≤-nonNeg (+ m) (ℤ.+≤+ (ℕₚ.n≤1+n (suc j)))))
