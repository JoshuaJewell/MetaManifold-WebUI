-- SPDX-License-Identifier: AGPL-3.0-only
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- The Benjamini–Hochberg step-up adjustment, over exact rationals, with the
-- properties a reader of a q-value table relies on proved rather than tested.
--
-- `bhAdjust` below models `bh_adjust` in `src/analysis/differential.jl`
-- (branch feat/differential-abundance-nb, commit 6c3e8b2):
--
--   Julia                                        here
--   ─────────────────────────────────────────    ──────────────────────────────
--   n = length(p)                                m = length ps
--   order = sortperm(p; rev=true)                order = sort (indexFrom 0 ps),
--                                                  ascending, lexicographic on
--                                                  (p, index)
--   running = 1.0                                running m j: the running
--   for (k, idx) in enumerate(order)               minimum, listed in ascending
--     rank = n - k + 1                             rank order; `minInf` starts
--     running = min(running, n/rank * p[idx])      it from the largest rank
--     out[idx] = min(running, 1.0)               stepUp = map (_⊓ 1ℚ) ∘ running
--                                                restore: output position idx
--                                                  reads the value whose key
--                                                  has index idx
--
-- Two differences.  The first is proved harmless; the second follows from a
-- proved lemma by a short argument about the Julia code, which is not itself
-- formalised:
--
--   * Julia seeds `running` with 1.0, so its running value at rank i is
--     min(1, min_{j ≥ i} n/j · p_(j)) and the final clamp is idempotent.
--     The model seeds with the first scaled value and clamps afterwards.
--     `stepUp-envelope` proves the model's values equal that same envelope,
--     so the two compute the same list.
--   * Tie-breaking.  Julia walks a stable DESCENDING sort, so of two equal
--     p-values the earlier input index gets the HIGHER rank; the model's
--     ascending (p, index) sort gives it the lower rank.  The sorted list of
--     p-values is the same either way, and tied p-values receive equal
--     adjusted values whatever their order (`bh-monotone` applied in both
--     directions), so which tied index reads which tied slot cannot change
--     the output.
--
-- Julia's ranks are one-based; here they are zero-based (`j = rank - 1`),
-- so Julia's `n / rank` is `rankFactor m j = m / (j + 1)` from the Prelude.
--
-- Julia rejects any p-value outside [0, 1] (and any non-finite value) with an
-- ArgumentError before it computes anything.  Over ℚ there is nothing
-- non-finite, so that guard becomes the hypothesis `All Valid ps` on every
-- theorem: the theorems speak about exactly the inputs Julia accepts.
--
-- Theorems (names as in docs/formal/agda-bh-scaling-proofs.md):
--
--   bh-length          the output has one entry per input
--   bh-nonneg          R-BH-1   every adjusted value is ≥ 0
--   bh-dominates-p     R-BH-2a  adjusted_i ≥ p_i
--   bh-at-most-one     R-BH-2b  adjusted_i ≤ 1
--   bh-monotone        R-BH-2c  p_i ≤ p_j ⇒ adjusted_i ≤ adjusted_j
--                              (ties included, any tie-break)
--   bh-envelope        R-BH-2d  the adjusted value at sorted rank k is
--                              min(1, min_{j ≥ k} p_(j) · m/(j+1))
--   bh-output-is-envelope       every (p, adjusted) pair of the output is a
--                              (p_(k), envelope_k) pair of the sorted stage
--
-- The archived residue file recorded "q_i ≥ p_i is FALSE in general".  That is
-- true of the bare scaling `p · M/(j+1)` for an arbitrary M (take M < j + 1),
-- but not of the procedure: the rank of a p-value never exceeds the family
-- size, so every factor m/(j+1) is ≥ 1, and R-BH-2a holds.  It is proved here.

{-# OPTIONS --safe --without-K #-}

module MetaManifold.BenjaminiHochberg where

open import MetaManifold.Prelude

open import Data.Integer.Base as ℤ using (ℤ; +_; +[1+_])
import Data.Integer.Properties as ℤₚ
open import Data.Nat.Base as ℕ using (ℕ; zero; suc; z≤n; s≤s)
import Data.Nat.Properties as ℕₚ
open import Data.List.Base using (List; []; _∷_; map; zip; length; foldr)
import Data.List.Properties as Listₚ
open import Data.List.Relation.Unary.All as All using (All; []; _∷_)
import Data.List.Relation.Unary.All.Properties as Allₚ
open import Data.List.Relation.Unary.Any as Any using (Any; here; there)
open import Data.List.Relation.Unary.Linked as Linked using (Linked; []; [-]; _∷_)
import Data.List.Relation.Unary.Linked.Properties as Linkedₚ
open import Data.List.Relation.Unary.AllPairs as AllPairs using (AllPairs; []; _∷_)
open import Data.List.Membership.Propositional using (_∈_)
import Data.List.Membership.Propositional.Properties as ∈ₚ
open import Data.List.Relation.Binary.Permutation.Propositional using (_↭_; ↭-sym)
import Data.List.Relation.Binary.Permutation.Propositional.Properties as ↭ₚ
open import Data.Product.Base using (_×_; _,_; proj₁; proj₂; ∃; ∃-syntax)
import Data.Product.Relation.Binary.Lex.NonStrict as Lex
open import Data.Sum.Base using (_⊎_; inj₁; inj₂)
open import Data.Empty using (⊥-elim)
open import Data.Rational.Base as ℚ
  using (ℚ; 0ℚ; 1ℚ; _*_; _⊓_; _≤_; toℚᵘ; fromℚᵘ)
import Data.Rational.Properties as ℚₚ
open import Data.Rational.Unnormalised.Base as ℚᵘ
  using (ℚᵘ; mkℚᵘ; ↥_; ↧ₙ_; *≡*)
  renaming (_≤_ to _≤ℚ_)
import Data.Rational.Unnormalised.Properties as ℚᵘₚ
open import Relation.Binary.Bundles using (DecTotalOrder)
open import Relation.Binary.PropositionalEquality
  using (_≡_; refl; sym; trans; cong; cong₂; subst)
open import Relation.Nullary using (yes; no)

------------------------------------------------------------------------
-- Part 1. The scaling step over ℚᵘ (ported from the archive)

-- `M/(j+1) · n/(d+1)`, written with the denominator expanded so that the
-- fraction is built by a single `mkℚᵘ` and carries no intermediate rounding.
-- `j` is a zero-based rank; `j + d + j * d` is `(j+1)(d+1) - 1`, the
-- denominator-minus-one that `mkℚᵘ` expects.
bhScale : ℕ → ℕ → ℤ → ℕ → ℚᵘ
bhScale M j n d = mkℚᵘ (+ M ℤ.* n) (j ℕ.+ d ℕ.+ j ℕ.* d)

-- The numerator of the scaled value is `M·n`.
bhScale-numerator : ∀ M j n d → ↥ (bhScale M j n d) ≡ + M ℤ.* n
bhScale-numerator M j n d = refl

-- The denominator of the scaled value is `(j+1)(d+1)`.
bhScale-denominator :
  ∀ M j n d → ↧ₙ (bhScale M j n d) ≡ suc (j ℕ.+ d ℕ.+ j ℕ.* d)
bhScale-denominator M j n d = refl

-- A bigger family can never produce a smaller scaled value.
monotone-in-family-size :
  ∀ M M′ j n d → M ℕ.≤ M′ → .{{_ : ℤ.NonNegative n}} →
  bhScale M j n d ≤ℚ bhScale M′ j n d
monotone-in-family-size M M′ j n d M≤M′ {{np}} =
  numerator-monotone-≤ (+ M ℤ.* n) (+ M′ ℤ.* n) (j ℕ.+ d ℕ.+ j ℕ.* d)
    (ℤₚ.≤-trans (ℤₚ.≤-reflexive (ℤₚ.*-comm (+ M) n))
      (ℤₚ.≤-trans (ℤₚ.*-monoˡ-≤-nonNeg n {{np}} (ℤ.+≤+ M≤M′))
                  (ℤₚ.≤-reflexive (ℤₚ.*-comm n (+ M′)))))

------------------------------------------------------------------------
-- Part 2. The procedure over ℚ

-- A p-value Julia accepts: a number in [0, 1].
Valid : ℚ → Set
Valid p = 0ℚ ≤ p × p ≤ 1ℚ

-- Julia's `n / rank * p[idx]` at zero-based rank j = rank - 1.
scaled : ℕ → ℕ → ℚ → ℚ
scaled m j p = p * rankFactor m j

-- Bridge to Part 1: the ℚ scaling is, after forgetting normalisation, exactly
-- the archived ℚᵘ `bhScale` applied to the p-value's numerator and
-- denominator-minus-one.  So Part 1's lemmas are lemmas about `scaled`.
scaled≃bhScale : ∀ m j p →
  toℚᵘ (scaled m j p) ℚᵘ.≃ bhScale m j (ℚ.numerator p) (ℚ.denominator-1 p)
scaled≃bhScale m j p@record{} = ℚᵘₚ.≃-trans (ℚₚ.toℚᵘ-homo-* p (rankFactor m j))
  (ℚᵘₚ.≃-trans (ℚᵘₚ.*-cong (ℚᵘₚ.≃-refl {toℚᵘ p}) (ℚₚ.toℚᵘ-fromℚᵘ (mkℚᵘ (+ m) j)))
    (*≡* (cong₂ ℤ._*_ (ℤₚ.*-comm (ℚ.numerator p) (+ m))
                      (cong +[1+_] (sym denom)))))
  where
  d = ℚ.denominator-1 p
  -- (d+1)(j+1) - 1 and (j+1)(d+1) - 1 agree.
  denom : j ℕ.+ d ℕ.* suc j ≡ j ℕ.+ d ℕ.+ j ℕ.* d
  denom = trans (cong (j ℕ.+_) (ℕₚ.*-suc d j))
           (trans (sym (ℕₚ.+-assoc j d (d ℕ.* j)))
                  (cong ((j ℕ.+ d) ℕ.+_) (ℕₚ.*-comm d j)))

-- The running minimum seeded with the first value folded in (the largest
-- rank).  Julia seeds with 1.0 instead; `stepUp-envelope` shows the two agree
-- once clamped.
minInf : ℚ → List ℚ → ℚ
minInf x []      = x
minInf x (r ∷ _) = r ⊓ x

-- The running minima of Julia's loop over ranks n, n-1, ..., listed in rank order, for a
-- sorted list whose first element sits at zero-based rank j.
running : ℕ → ℕ → List ℚ → List ℚ
running m j []       = []
running m j (p ∷ ps) =
  minInf (scaled m j p) (running m (suc j) ps) ∷ running m (suc j) ps

-- Julia's `min(running, 1.0)`.
stepUp : ℕ → List ℚ → List ℚ
stepUp m ps = map (_⊓ 1ℚ) (running m 0 ps)

-- Pair each value with its (zero-based) position: Julia's index `idx`.
indexFrom : ℕ → List ℚ → List (ℚ × ℕ)
indexFrom k []       = []
indexFrom k (p ∷ ps) = (p , k) ∷ indexFrom (suc k) ps

-- Ascending order on (p, index).  Julia sorts descending (stably); the
-- tie order differs and is harmless (see the header).
keyOrder : DecTotalOrder _ _ _
keyOrder = Lex.×-decTotalOrder ℚₚ.≤-decTotalOrder ℕₚ.≤-decTotalOrder

open import Data.List.Sort.MergeSort keyOrder using (sort; sort-↭; sort-↗)

-- `out[idx] = min(running, 1.0)` read from the other side: the value at output
-- position i is the adjusted value of the first sorted
-- entry whose index is i.  The `0ℚ` default is never reached (`restore-spec`),
-- as Julia's `undef` slots are never left unwritten.
restore : List ((ℚ × ℕ) × ℚ) → ℕ → ℚ
restore []                    i = 0ℚ
restore (((p , k) , a) ∷ ts) i with k ℕₚ.≟ i
... | yes _ = a
... | no  _ = restore ts i

-- The intermediate values of `bh_adjust`, named as in Julia.
module Pipeline (ps : List ℚ) where
  m : ℕ
  m = length ps

  keyed : List (ℚ × ℕ)
  keyed = indexFrom 0 ps

  order : List (ℚ × ℕ)
  order = sort keyed

  sortedP : List ℚ
  sortedP = map proj₁ order

  adjusted : List ℚ
  adjusted = stepUp m sortedP

  triples : List ((ℚ × ℕ) × ℚ)
  triples = zip order adjusted

-- The Benjamini–Hochberg adjusted p-values, in input order.
bhAdjust : List ℚ → List ℚ
bhAdjust ps = map (λ k → restore triples (proj₂ k)) keyed
  where open Pipeline ps

------------------------------------------------------------------------
-- Part 3. The sorted stage: the envelope form (R-BH-2d)

-- The scaled values p_(j) · m/(j+1) of a sorted list starting at rank j.
scaledFrom : ℕ → ℕ → List ℚ → List ℚ
scaledFrom m j []       = []
scaledFrom m j (p ∷ ps) = scaled m j p ∷ scaledFrom m (suc j) ps

-- For each suffix, min(1, min of the suffix).
envelopes : List ℚ → List ℚ
envelopes []       = []
envelopes (x ∷ xs) = foldr _⊓_ 1ℚ (x ∷ xs) ∷ envelopes xs

-- The head of the clamped running minimum is the clamped suffix minimum.
private
  head-eq : ∀ m j p ps →
    minInf (scaled m j p) (running m (suc j) ps) ⊓ 1ℚ
      ≡ scaled m j p ⊓ foldr _⊓_ 1ℚ (scaledFrom m (suc j) ps)
  head-eq m j p []       = refl
  head-eq m j p (q ∷ qs) =
    trans (cong (_⊓ 1ℚ) (ℚₚ.⊓-comm r s))
      (trans (ℚₚ.⊓-assoc s r 1ℚ)
             (cong (s ⊓_) (head-eq m (suc j) q qs)))
    where
    s = scaled m j p
    r = minInf (scaled m (suc j) q) (running m (suc (suc j)) qs)

-- R-BH-2d, sorted-stage form: the Julia loop computes the envelope.
stepUp-envelope : ∀ m j ps →
  map (_⊓ 1ℚ) (running m j ps) ≡ envelopes (scaledFrom m j ps)
stepUp-envelope m j []       = refl
stepUp-envelope m j (p ∷ ps) =
  cong₂ _∷_ (head-eq m j p ps) (stepUp-envelope m (suc j) ps)

-- R-BH-2d for the pipeline: the adjusted values at sorted rank k are
-- min(1, min_{j ≥ k} p_(j) · m/(j+1)).
bh-envelope : ∀ ps →
  Pipeline.adjusted ps
    ≡ envelopes (scaledFrom (length ps) 0 (Pipeline.sortedP ps))
bh-envelope ps = stepUp-envelope (length ps) 0 (Pipeline.sortedP ps)

------------------------------------------------------------------------
-- Part 4. Facts about the scaling and the fold

private
  0≤1 : 0ℚ ≤ 1ℚ
  0≤1 = ℚₚ.toℚᵘ-cancel-≤ (ℚᵘ.*≤* (ℤ.+≤+ z≤n))

  -- Scaling inside the family never lowers a non-negative p-value.
  scaled-≥ : ∀ m j p → suc j ℕ.≤ m → 0ℚ ≤ p → p ≤ scaled m j p
  scaled-≥ m j p j<m 0≤p =
    ℚₚ.≤-trans (ℚₚ.≤-reflexive (sym (ℚₚ.*-identityʳ p)))
      (ℚₚ.*-monoˡ-≤-nonNeg p {{ℚ.nonNegative 0≤p}} (rankFactor-≥1 m j j<m))

  -- Scaling is monotone in the p-value.
  scaled-mono : ∀ m j {p q} → p ≤ q → scaled m j p ≤ scaled m j q
  scaled-mono m j p≤q =
    ℚₚ.*-monoʳ-≤-nonNeg (rankFactor m j)
      {{ℚ.nonNegative (rankFactor-nonNeg m j)}} p≤q

  -- The same non-negative p-value scaled at a later rank is no larger.
  scaled-anti : ∀ m j p → 0ℚ ≤ p → scaled m (suc j) p ≤ scaled m j p
  scaled-anti m j p 0≤p =
    ℚₚ.*-monoˡ-≤-nonNeg p {{ℚ.nonNegative 0≤p}} (rankFactor-anti m j)

  -- A clamped fold never exceeds one.
  fold-≤1 : ∀ xs → foldr _⊓_ 1ℚ xs ≤ 1ℚ
  fold-≤1 []       = ℚₚ.≤-refl
  fold-≤1 (x ∷ xs) = ℚₚ.≤-trans (ℚₚ.p⊓q≤q x _) (fold-≤1 xs)

  -- A lower bound of one and of every scaled value bounds the fold below.
  fold-lower : ∀ m j ps q → q ≤ 1ℚ → All (q ≤_) ps → All (0ℚ ≤_) ps →
    j ℕ.+ length ps ℕ.≤ m → q ≤ foldr _⊓_ 1ℚ (scaledFrom m j ps)
  fold-lower m j []       q q≤1 []          []          _ = q≤1
  fold-lower m j (p ∷ ps) q q≤1 (q≤p ∷ q≤s) (0≤p ∷ 0≤s) bound =
    ℚₚ.⊓-glb (ℚₚ.≤-trans q≤p (scaled-≥ m j p j<m 0≤p))
             (fold-lower m (suc j) ps q q≤1 q≤s 0≤s bound′)
    where
    bound′ : suc j ℕ.+ length ps ℕ.≤ m
    bound′ = ℕₚ.≤-trans (ℕₚ.≤-reflexive (sym (ℕₚ.+-suc j (length ps)))) bound
    j<m : suc j ℕ.≤ m
    j<m = ℕₚ.≤-trans (ℕₚ.m≤m+n (suc j) (length ps)) bound′

------------------------------------------------------------------------
-- Part 5. The sorted stage: bounds (R-BH-1, 2a, 2b)

-- The sorted stage: each sorted p-value next to its adjusted value.
Z : ℕ → ℕ → List ℚ → List (ℚ × ℚ)
Z m j ps = zip ps (envelopes (scaledFrom m j ps))

-- The three bounds, as one predicate on (p, adjusted).
Bounded : ℚ × ℚ → Set
Bounded (p , a) = 0ℚ ≤ a × p ≤ a × a ≤ 1ℚ

private
  -- In a sorted list, the head is below everything.
  head-below : ∀ {p ps} → Linked _≤_ (p ∷ ps) → All (p ≤_) (p ∷ ps)
  head-below lnk = Linkedₚ.Linked⇒All ℚₚ.≤-trans ℚₚ.≤-refl lnk

  Z-bounded : ∀ m j ps → Linked _≤_ ps → All Valid ps →
    j ℕ.+ length ps ℕ.≤ m → All Bounded (Z m j ps)
  Z-bounded m j []       _   []                _     = []
  Z-bounded m j (p ∷ ps) lnk (vp@(0≤p , p≤1) ∷ vs) bound =
    ( fold-lower m j (p ∷ ps) 0ℚ 0≤1 nonneg nonneg bound
    , fold-lower m j (p ∷ ps) p p≤1 (head-below lnk) nonneg bound
    , fold-≤1 (scaledFrom m j (p ∷ ps)) )
    ∷ Z-bounded m (suc j) ps (Linked.tail lnk) vs bound′
    where
    nonneg : All (0ℚ ≤_) (p ∷ ps)
    nonneg = All.map proj₁ (vp ∷ vs)
    bound′ : suc j ℕ.+ length ps ℕ.≤ m
    bound′ = ℕₚ.≤-trans (ℕₚ.≤-reflexive (sym (ℕₚ.+-suc j (length ps)))) bound

------------------------------------------------------------------------
-- Part 6. The sorted stage: monotonicity, ties included (R-BH-2c)

-- Along the sorted stage: p never decreases, the adjusted value never
-- decreases, and where p ties the adjusted value ties too.
Step : ℚ × ℚ → ℚ × ℚ → Set
Step (p , a) (p′ , a′) = p ≤ p′ × a ≤ a′ × (p′ ≤ p → a′ ≤ a)

private
  Step-trans : ∀ {x y z} → Step x y → Step y z → Step x z
  Step-trans (p≤q , a≤b , tie₁) (q≤r , b≤c , tie₂) =
    ( ℚₚ.≤-trans p≤q q≤r
    , ℚₚ.≤-trans a≤b b≤c
    , λ r≤p → ℚₚ.≤-trans (tie₂ (ℚₚ.≤-trans r≤p p≤q))
                         (tie₁ (ℚₚ.≤-trans q≤r r≤p)) )

  -- Adjacent sorted entries are related by Step.  The tie case is where
  -- the cumulative minimum is needed: without it the later entry's
  -- envelope need not lie below the earlier entry's scaled value.
  Z-linked : ∀ m j ps → Linked _≤_ ps → All (0ℚ ≤_) ps →
    Linked Step (Z m j ps)
  Z-linked m j []            _           _                = []
  Z-linked m j (p ∷ [])      _           _                = [-]
  Z-linked m j (p ∷ q ∷ ps) (p≤q ∷ lnk) (0≤p ∷ 0≤q ∷ nn) =
    ( p≤q
    , ℚₚ.p⊓q≤q s F
    , (λ q≤p → ℚₚ.⊓-glb (F≤s q≤p) ℚₚ.≤-refl) )
    ∷ Z-linked m (suc j) (q ∷ ps) lnk (0≤q ∷ nn)
    where
    s = scaled m j p
    t = scaled m (suc j) q
    F = foldr _⊓_ 1ℚ (scaledFrom m (suc j) (q ∷ ps))
    F≤s : q ≤ p → F ≤ s
    F≤s q≤p = ℚₚ.≤-trans (ℚₚ.p⊓q≤p t _)
                (ℚₚ.≤-trans (scaled-mono m (suc j) q≤p) (scaled-anti m j p 0≤p))

  -- Any two members of an all-pairs-related list are equal or related.
  pairwise : ∀ {R : ℚ × ℚ → ℚ × ℚ → Set} {xs x y} → AllPairs R xs →
    x ∈ xs → y ∈ xs → x ≡ y ⊎ (R x y ⊎ R y x)
  pairwise (rx ∷ rs) (here refl) (here refl) = inj₁ refl
  pairwise (rx ∷ rs) (here refl) (there y∈)  = inj₂ (inj₁ (All.lookup rx y∈))
  pairwise (rx ∷ rs) (there x∈)  (here refl) = inj₂ (inj₂ (All.lookup rx x∈))
  pairwise (rx ∷ rs) (there x∈)  (there y∈)  = pairwise rs x∈ y∈

  Z-monotone : ∀ m j ps → Linked _≤_ ps → All (0ℚ ≤_) ps →
    ∀ {x y} → x ∈ Z m j ps → y ∈ Z m j ps →
    proj₁ x ≤ proj₁ y → proj₂ x ≤ proj₂ y
  Z-monotone m j ps lnk nn x∈ y∈ px≤py
    with pairwise (Linkedₚ.Linked⇒AllPairs Step-trans (Z-linked m j ps lnk nn)) x∈ y∈
  ... | inj₁ refl                   = ℚₚ.≤-refl
  ... | inj₂ (inj₁ (_ , a≤b , _))   = a≤b
  ... | inj₂ (inj₂ (_ , _ , tie))   = tie px≤py

------------------------------------------------------------------------
-- Part 7. From the sorted stage back to input order

private
  map-proj₁-indexFrom : ∀ k ps → map proj₁ (indexFrom k ps) ≡ ps
  map-proj₁-indexFrom k []       = refl
  map-proj₁-indexFrom k (p ∷ ps) = cong (p ∷_) (map-proj₁-indexFrom (suc k) ps)

  length-indexFrom : ∀ k ps → length (indexFrom k ps) ≡ length ps
  length-indexFrom k []       = refl
  length-indexFrom k (p ∷ ps) = cong suc (length-indexFrom (suc k) ps)

  All-indexFrom : ∀ {P : ℚ → Set} k ps → All P ps →
    All (λ x → P (proj₁ x)) (indexFrom k ps)
  All-indexFrom k []       []         = []
  All-indexFrom k (p ∷ ps) (pp ∷ pps) = pp ∷ All-indexFrom (suc k) ps pps

  indexFrom-≥ : ∀ {a i} k ps → (a , i) ∈ indexFrom k ps → k ℕ.≤ i
  indexFrom-≥ k (p ∷ ps) (here refl) = ℕₚ.≤-refl
  indexFrom-≥ k (p ∷ ps) (there x∈)  = ℕₚ.≤-trans (ℕₚ.n≤1+n k) (indexFrom-≥ (suc k) ps x∈)

  -- An index determines its p-value.
  indexFrom-functional : ∀ {a b i} k ps →
    (a , i) ∈ indexFrom k ps → (b , i) ∈ indexFrom k ps → a ≡ b
  indexFrom-functional k (p ∷ ps) (here a≡) (here b≡) =
    trans (cong proj₁ a≡) (sym (cong proj₁ b≡))
  indexFrom-functional k (p ∷ ps) (here a≡) (there b∈)  =
    ⊥-elim (ℕₚ.<-irrefl (sym (cong proj₂ a≡)) (indexFrom-≥ (suc k) ps b∈))
  indexFrom-functional k (p ∷ ps) (there a∈)  (here b≡) =
    ⊥-elim (ℕₚ.<-irrefl (sym (cong proj₂ b≡)) (indexFrom-≥ (suc k) ps a∈))
  indexFrom-functional k (p ∷ ps) (there a∈)  (there b∈)  =
    indexFrom-functional (suc k) ps a∈ b∈

  length-running : ∀ m j ps → length (running m j ps) ≡ length ps
  length-running m j []       = refl
  length-running m j (p ∷ ps) = cong suc (length-running m (suc j) ps)

  map-proj₁-zip : ∀ {A B : Set} (xs : List A) (ys : List B) →
    length xs ≡ length ys → map proj₁ (zip xs ys) ≡ xs
  map-proj₁-zip []       []       _  = refl
  map-proj₁-zip (x ∷ xs) (y ∷ ys) eq =
    cong (x ∷_) (map-proj₁-zip xs ys (ℕₚ.suc-injective eq))

  zip-map-map : ∀ {A B C : Set} (f : A → B) (g : A → C) xs →
    zip (map f xs) (map g xs) ≡ map (λ x → f x , g x) xs
  zip-map-map f g []       = refl
  zip-map-map f g (x ∷ xs) = cong ((f x , g x) ∷_) (zip-map-map f g xs)

  zip-map-left : ∀ {A B C : Set} (f : A → B) (xs : List A) (ys : List C) →
    zip (map f xs) ys ≡ map (λ t → f (proj₁ t) , proj₂ t) (zip xs ys)
  zip-map-left f []       ys       = refl
  zip-map-left f (x ∷ xs) []       = refl
  zip-map-left f (x ∷ xs) (y ∷ ys) = cong ((f x , y) ∷_) (zip-map-left f xs ys)

  -- `restore` finds the adjusted value of an entry with the requested index.
  restore-spec : ∀ ts i → Any (λ t → proj₂ (proj₁ t) ≡ i) ts →
    ∃[ t ] (t ∈ ts × proj₂ (proj₁ t) ≡ i × restore ts i ≡ proj₂ t)
  restore-spec (((p , k) , a) ∷ ts) i any with k ℕₚ.≟ i
  ... | yes k≡i = ((p , k) , a) , here refl , k≡i , refl
  restore-spec (((p , k) , a) ∷ ts) i (here k≡i) | no k≢i = ⊥-elim (k≢i k≡i)
  restore-spec (((p , k) , a) ∷ ts) i (there any) | no k≢i
    with restore-spec ts i any
  ... | t , t∈ , idx , eq = t , there t∈ , idx , eq

module _ (ps : List ℚ) where
  open Pipeline ps

  -- The sort neither drops nor duplicates entries.
  length-sortedP : length sortedP ≡ m
  length-sortedP =
    trans (Listₚ.length-map proj₁ order)
      (trans (↭ₚ.↭-length (sort-↭ keyed)) (length-indexFrom 0 ps))

  private
    length-adjusted : length adjusted ≡ length order
    length-adjusted =
      trans (Listₚ.length-map (_⊓ 1ℚ) (running m 0 sortedP))
        (trans (length-running m 0 sortedP) (Listₚ.length-map proj₁ order))

    map-proj₁-triples : map proj₁ triples ≡ order
    map-proj₁-triples = map-proj₁-zip order adjusted (sym length-adjusted)

    order↭keyed : order ↭ keyed
    order↭keyed = sort-↭ keyed

    -- Lexicographic order on (p, index) implies order on p.
    lex⇒≤ : ∀ {x y} → DecTotalOrder._≤_ keyOrder x y → proj₁ x ≤ proj₁ y
    lex⇒≤ (inj₁ (p≤q , _)) = p≤q
    lex⇒≤ (inj₂ (p≡q , _)) = ℚₚ.≤-reflexive p≡q

  -- The sorted p-values really are sorted.
  sortedP-sorted : Linked _≤_ sortedP
  sortedP-sorted = Linkedₚ.map⁺ (Linked.map lex⇒≤ (sort-↗ keyed))

  -- Validity survives the sort.
  sortedP-valid : All Valid ps → All Valid sortedP
  sortedP-valid vs = Allₚ.map⁺
    (↭ₚ.All-resp-↭ (↭-sym order↭keyed) (All-indexFrom 0 ps vs))

  -- The sorted stage of the pipeline is `Z m 0 sortedP`.
  zip-sorted-adjusted : zip sortedP adjusted ≡ Z m 0 sortedP
  zip-sorted-adjusted = cong (zip sortedP) (bh-envelope ps)

  -- Every (p, adjusted) pair of the output is a (p_(k), envelope_k) pair of
  -- the sorted stage.  This is the whole content of Julia's scatter-back
  -- `out[idx] = ...` that the theorems need.
  bh-output-is-envelope : ∀ {x} → x ∈ zip ps (bhAdjust ps) → x ∈ Z m 0 sortedP
  bh-output-is-envelope {x} x∈ = subst (x ∈_) zip-sorted-adjusted x∈Zs
    where
    g : ℚ × ℕ → ℚ
    g k = restore triples (proj₂ k)

    out-shape : zip ps (bhAdjust ps) ≡ map (λ k → proj₁ k , g k) keyed
    out-shape = trans
      (cong (λ l → zip l (bhAdjust ps)) (sym (map-proj₁-indexFrom 0 ps)))
      (zip-map-map proj₁ g keyed)

    x∈Zs : x ∈ zip sortedP adjusted
    x∈Zs with ∈ₚ.∈-map⁻ (λ k → proj₁ k , g k) (subst (x ∈_) out-shape x∈)
    ... | k , k∈keyed , refl
      with ∈ₚ.∈-map⁻ proj₁
             (subst (k ∈_) (sym map-proj₁-triples)
               (↭ₚ.∈-resp-↭ (↭-sym order↭keyed) k∈keyed))
    ... | t₀ , t₀∈ , refl
      with restore-spec triples (proj₂ (proj₁ t₀))
             (Any.map (λ e → cong (λ u → proj₂ (proj₁ u)) (sym e)) t₀∈)
    ... | t , t∈ , idx , eq = subst (_∈ zip sortedP adjusted) pair-eq t-in
      where
      key∈keyed : ∀ {u} → u ∈ triples → proj₁ u ∈ keyed
      key∈keyed u∈ = ↭ₚ.∈-resp-↭ order↭keyed
        (subst (_ ∈_) map-proj₁-triples (∈ₚ.∈-map⁺ proj₁ u∈))

      same-p : proj₁ (proj₁ t) ≡ proj₁ (proj₁ t₀)
      same-p = indexFrom-functional 0 ps
        (subst (λ i → (proj₁ (proj₁ t) , i) ∈ keyed) idx (key∈keyed t∈))
        (key∈keyed t₀∈)

      t-in : (proj₁ (proj₁ t) , proj₂ t) ∈ zip sortedP adjusted
      t-in = subst ((proj₁ (proj₁ t) , proj₂ t) ∈_)
               (sym (zip-map-left proj₁ order adjusted))
               (∈ₚ.∈-map⁺ (λ u → proj₁ (proj₁ u) , proj₂ u) t∈)

      pair-eq : (proj₁ (proj₁ t) , proj₂ t)
                ≡ (proj₁ (proj₁ t₀) , restore triples (proj₂ (proj₁ t₀)))
      pair-eq = cong₂ _,_ same-p (sym eq)

------------------------------------------------------------------------
-- Part 8. The theorems, in input order

-- The output has one entry per input.
bh-length : ∀ ps → length (bhAdjust ps) ≡ length ps
bh-length ps = trans (Listₚ.length-map _ (indexFrom 0 ps)) (length-indexFrom 0 ps)

-- All three bounds at once, for every (p_i, adjusted_i).
bh-bounded : ∀ ps → All Valid ps → All Bounded (zip ps (bhAdjust ps))
bh-bounded ps vs = All.tabulate λ x∈ →
  All.lookup
    (Z-bounded m 0 sortedP (sortedP-sorted ps) (sortedP-valid ps vs)
      (ℕₚ.≤-reflexive (length-sortedP ps)))
    (bh-output-is-envelope ps x∈)
  where open Pipeline ps

-- R-BH-1: every adjusted value is non-negative.
bh-nonneg : ∀ ps → All Valid ps →
  All (λ x → 0ℚ ≤ proj₂ x) (zip ps (bhAdjust ps))
bh-nonneg ps vs = All.map proj₁ (bh-bounded ps vs)

-- R-BH-2a: every adjusted value is at least its own p-value.
bh-dominates-p : ∀ ps → All Valid ps →
  All (λ x → proj₁ x ≤ proj₂ x) (zip ps (bhAdjust ps))
bh-dominates-p ps vs = All.map (λ b → proj₁ (proj₂ b)) (bh-bounded ps vs)

-- R-BH-2b: every adjusted value is at most one.
bh-at-most-one : ∀ ps → All Valid ps →
  All (λ x → proj₂ x ≤ 1ℚ) (zip ps (bhAdjust ps))
bh-at-most-one ps vs = All.map (λ b → proj₂ (proj₂ b)) (bh-bounded ps vs)

-- R-BH-2c: a smaller (or equal) p-value never gets a larger adjusted value.
-- Stated for any two input positions, so ties (equal p-values at different
-- positions) are covered in both directions: they get equal adjusted values,
-- whatever order the sort put them in.
bh-monotone : ∀ ps → All Valid ps →
  ∀ {x y} → x ∈ zip ps (bhAdjust ps) → y ∈ zip ps (bhAdjust ps) →
  proj₁ x ≤ proj₁ y → proj₂ x ≤ proj₂ y
bh-monotone ps vs x∈ y∈ =
  Z-monotone m 0 sortedP (sortedP-sorted ps)
    (All.map proj₁ (sortedP-valid ps vs))
    (bh-output-is-envelope ps x∈) (bh-output-is-envelope ps y∈)
  where open Pipeline ps
