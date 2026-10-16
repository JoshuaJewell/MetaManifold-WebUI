-- SPDX-License-Identifier: AGPL-3.0-only
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- PERMANOVA's permutation bookkeeping: the two things `run_permanova` in
-- `src/analysis/analysis.jl` decides for itself rather than taking from vegan
-- (branch fix/permanova-terms, PR #35), proved.
--
--   1. `min_p_value = 1 / (n_perm + 1)`.  vegan's permutation p-value is
--      (k + 1)/(B + 1), where B is the number of permutations performed and k
--      how many of them gave a statistic at least the observed one (0 ≤ k ≤ B).
--      The theorems say that 1/(B + 1) is the smallest p-value B permutations
--      can produce, that it is attained, and that no threshold below it can be
--      reached.
--   2. Which terms have their p-value withheld.  With permutations restricted
--      to blocks, Julia withholds the p-value of a term when
--      `_varies_within_any_block(values, blocks)` is false.  The theorems say
--      that this happens exactly when every within-block permutation leaves the
--      term's labels where they were, so the permutations never exchange them.
--      They also say that when it is true, some within-block permutation does
--      exchange them. This characterises the block-constancy filter only;
--      it does not establish statistical testability or a reported p-value.
--
-- What is NOT modelled, and why:
--
--   * The pseudo-F statistic and vegan's distribution of it.  The p-value
--     theorems take k as given; they are about the arithmetic of (k + 1)/(B + 1),
--     not about which permutations exceed the observed F.
--   * How many within-block permutations exist.  B is whatever `permute`
--     reports (`attr(perm_res, "control")$nperm`).  Six blocks of two samples
--     admit 2⁶ = 64 arrangements, the identity and 63 others, which is the 63 in
--     test_permanova.jl.  That count is the product of the blocks' factorials;
--     it is not proved here.
--   * Unblocked runs.  Julia only asks `_varies_within_any_block` when some
--     block label repeats (`blocked`).  Without repeated blocks the
--     permutations are unrestricted and nothing is withheld; that branch is
--     Julia's `if`, not a theorem.
--   * Julia's `isequal` on NaN or `missing`.  Equality of labels and of blocks
--     is an arbitrary decidable equality here.
--   * Floating point.  `1 / (n_perm + 1)` is exact here and a Float64 in Julia.
--
-- Theorems (names as in proofs/agda/README.md):
--
--   minP-attained          (0 + 1)/(B + 1) is 1/(B + 1)
--   p-floor                every permutation p-value is at least 1/(B + 1)
--   p-at-most-one          ... and at most 1 when k ≤ B
--   below-floor-unreachable  a threshold below 1/(B + 1) is below every p-value
--   minP-anti              more permutations never raise the floor
--   minP-63                with 63 permutations the floor is 1/64
--   no-p-at-1%-with-63     with 63 permutations no p-value is at or below 0.01
--   never-exchanged        a block-constant term is fixed by every
--                          within-block permutation
--   exchanged              a term with two labels in one block is moved by a
--                          within-block permutation (a transposition)
--   untestable⇔fixed       block-constant ⇔ fixed by every within-block
--                          permutation
--   varies-false⇔Constant  the model of `_varies_within_any_block` returns
--                          false exactly on block-constant input
--   withheld⇔never-exchanged  ... so its false answer, which withholds the
--                          p-value, holds exactly when no within-block
--                          permutation exchanges the term's labels
--   tested⇔exchanged       ... and its true answer, which passes this filter,
--                          holds exactly when some within-block permutation
--                          does
--   constant-in-blocks, varies-in-a-block
--                          the two direct cases in test_permanova.jl
--   fixture-group-withheld, fixture-run-tested
--                          the twelve-sample fixture in test_permanova.jl: group
--                          is never exchanged, run is

{-# OPTIONS --safe --without-K #-}

module MetaManifold.Permanova where

open import Data.Bool.Base using (Bool; true; false)
open import Data.Empty using (⊥-elim)
open import Data.Fin.Base using (Fin)
open import Data.Fin.Properties using (_≟_)
open import Data.Fin.Permutation using (Permutation′; _⟨$⟩ʳ_; transpose)
open import Data.Integer.Base as ℤ using (+_; +[1+_])
import Data.Integer.Properties as ℤₚ
open import Data.List.Base using (List; []; _∷_; _++_; [_]; length; lookup)
open import Data.List.Membership.Propositional using (_∈_)
open import Data.List.Membership.Propositional.Properties
  using (∈-lookup; ∈-++⁺ˡ; ∈-++⁺ʳ; ∈-++⁻)
open import Data.List.Properties using (++-assoc; ++-identityʳ)
open import Data.List.Relation.Unary.Any using (here; index)
open import Data.List.Relation.Unary.Any.Properties using (lookup-index)
open import Data.Maybe.Base using (Maybe; just; nothing)
open import Data.Nat.Base as ℕ using (ℕ; suc; z≤n; s≤s)
import Data.Nat.Properties as ℕₚ
open import Data.Product.Base using (Σ; ∃₂; _×_; _,_; proj₁; proj₂)
open import Data.Rational.Base as ℚ using (ℚ; _/_; 1ℚ; _≤_; _<_)
import Data.Rational.Properties as ℚₚ
open import Data.Rational.Unnormalised.Base using (mkℚᵘ; *≤*)
open import Data.Sum.Base using (inj₁; inj₂)
open import Function.Bundles using (_⇔_; mk⇔; Equivalence)
open import Relation.Binary.Definitions using (DecidableEquality)
open import Relation.Binary.PropositionalEquality
  using (_≡_; _≢_; refl; sym; trans; cong; subst)
open import Relation.Nullary using (¬_; yes; no)
open import Relation.Nullary.Decidable using (toWitness)
open import Relation.Nullary.Negation using (contradiction)

open import MetaManifold.Prelude using (numerator-monotone-≤; fromℚᵘ-mono-≤)

------------------------------------------------------------------------
-- 1. The smallest p-value B permutations can produce

-- vegan's permutation p-value after B permutations, k of which gave a
-- statistic at least the observed one.
permP : ℕ → ℕ → ℚ
permP B k = + suc k / suc B

-- Julia's `min_p_value = 1 / (Int(n_perm) + 1)`.
minP : ℕ → ℚ
minP B = + 1 / suc B

-- The floor is a p-value: the one reported when no permutation reaches the
-- observed statistic.
minP-attained : ∀ B → permP B 0 ≡ minP B
minP-attained B = refl

-- No permutation p-value is below the floor.
p-floor : ∀ B k → minP B ≤ permP B k
p-floor B k = fromℚᵘ-mono-≤ {mkℚᵘ (+ 1) B} {mkℚᵘ (+ suc k) B}
  (numerator-monotone-≤ (+ 1) (+ suc k) B (ℤ.+≤+ (s≤s z≤n)))

-- With k ≤ B, none is above 1.
p-at-most-one : ∀ B k → k ℕ.≤ B → permP B k ≤ 1ℚ
p-at-most-one B k k≤B = fromℚᵘ-mono-≤ {mkℚᵘ (+ suc k) B} {mkℚᵘ (+ 1) 0}
  (*≤* (ℤₚ.≤-trans (ℤₚ.≤-reflexive (ℤₚ.*-identityʳ (+ suc k)))
         (ℤₚ.≤-trans (ℤ.+≤+ (s≤s k≤B))
                     (ℤₚ.≤-reflexive (sym (ℤₚ.*-identityˡ +[1+ B ]))))))

-- A significance threshold below the floor is below every p-value, so no term
-- can be significant at it, whatever the data.
below-floor-unreachable : ∀ B k {α} → α < minP B → α < permP B k
below-floor-unreachable B k α<floor = ℚₚ.<-≤-trans α<floor (p-floor B k)

-- More permutations never raise the floor.
minP-anti : ∀ {B B′} → B ℕ.≤ B′ → minP B′ ≤ minP B
minP-anti {B} {B′} B≤B′ = fromℚᵘ-mono-≤ {mkℚᵘ (+ 1) B′} {mkℚᵘ (+ 1) B}
  (*≤* (ℤₚ.*-monoˡ-≤-nonNeg (+ 1) (ℤ.+≤+ (s≤s B≤B′))))

-- The case test_permanova.jl measures: six individuals of two samples give 63
-- within-block permutations, so the floor is 1/64 and a 1% threshold can never
-- be met.
minP-63 : minP 63 ≡ + 1 / 64
minP-63 = refl

no-p-at-1%-with-63 : ∀ k → + 1 / 100 < permP 63 k
no-p-at-1%-with-63 k =
  below-floor-unreachable 63 k (toWitness {a? = (+ 1 / 100) ℚₚ.<? minP 63} _)

------------------------------------------------------------------------
-- 2. Within-block permutations and a term's labels

-- `n` samples; `block i` is sample i's block (the individual it came from) and
-- `label i` its value for one term of the model.
module Labels {A B : Set} {n : ℕ} (block : Fin n → B) (label : Fin n → A) where

  -- π keeps every sample in its own block: the only permutations `permute`
  -- draws when it is given `how(blocks = ...)`.
  WithinBlock : Permutation′ n → Set
  WithinBlock π = ∀ k → block (π ⟨$⟩ʳ k) ≡ block k

  -- π leaves the term's labels where they were.
  Fixes : Permutation′ n → Set
  Fixes π = ∀ k → label (π ⟨$⟩ʳ k) ≡ label k

  -- The term takes one value inside each block.
  BlockConstant : Set
  BlockConstant = ∀ i j → block i ≡ block j → label i ≡ label j

  -- A block-constant term's labels are never exchanged: every within-block
  -- permutation fixes them, so the permutation distribution is not a null
  -- distribution for that term.
  never-exchanged : BlockConstant → ∀ π → WithinBlock π → Fixes π
  never-exchanged bc π wb k = bc (π ⟨$⟩ʳ k) k (wb k)

  -- Swapping two samples of one block is a within-block permutation.
  transpose-within : ∀ {i j} → block i ≡ block j → WithinBlock (transpose i j)
  transpose-within {i} {j} bi≡bj k with k ≟ i
  ... | yes refl = sym bi≡bj
  ... | no _ with k ≟ j
  ...   | yes refl = bi≡bj
  ...   | no _ = refl

  transpose-at-i : ∀ (i j : Fin n) → transpose i j ⟨$⟩ʳ i ≡ j
  transpose-at-i i j with i ≟ i
  ... | yes _ = refl
  ... | no i≢i = contradiction refl i≢i

  -- Two different labels inside one block are exchanged by swapping them.
  exchanged : ∀ {i j} → block i ≡ block j → label i ≢ label j →
              Σ (Permutation′ n) λ π → WithinBlock π × ¬ Fixes π
  exchanged {i} {j} bi≡bj li≢lj =
    transpose i j , transpose-within bi≡bj ,
    λ fixes → li≢lj (sym (trans (sym (cong label (transpose-at-i i j))) (fixes i)))

  -- So "constant within every block" is exactly "never exchanged".
  untestable⇔fixed : BlockConstant ⇔ (∀ π → WithinBlock π → Fixes π)
  untestable⇔fixed = mk⇔ never-exchanged from
    where
    from : (∀ π → WithinBlock π → Fixes π) → BlockConstant
    from all-fixed i j bi≡bj =
      sym (trans (sym (cong label (transpose-at-i i j)))
                 (all-fixed (transpose i j) (transpose-within bi≡bj) i))

------------------------------------------------------------------------
-- 3. `_varies_within_any_block`, modelled and proved
--
-- Julia:
--
--   first_seen = Dict{eltype(blocks),eltype(values)}()
--   for (v, b) in zip(values, blocks)
--       isequal(get!(first_seen, b, v), v) || return true
--   end
--   false
--
-- A list of (value, block) pairs stands in for `zip(values, blocks)`, and an
-- association list for the dictionary.  `get!` returns the stored value, or
-- stores `v` and returns it; storing happens only on the `nothing` branch, as
-- in Julia, where the comparison that follows is `isequal(v, v)`.

module Varies {A B : Set} (_≟ᴬ_ : DecidableEquality A) (_≟ᴮ_ : DecidableEquality B) where

  -- `get(first_seen, b, nothing)`.
  firstSeen : List (B × A) → B → Maybe A
  firstSeen [] b = nothing
  firstSeen ((b′ , v) ∷ seen) b with b ≟ᴮ b′
  ... | yes _ = just v
  ... | no _ = firstSeen seen b

  go : List (B × A) → List (A × B) → Bool
  go seen [] = false
  go seen ((v , b) ∷ xs) with firstSeen seen b
  ... | nothing = go ((b , v) ∷ seen) xs
  ... | just w with w ≟ᴬ v
  ...   | yes _ = go seen xs
  ...   | no _ = true

  variesWithinAnyBlock : List (A × B) → Bool
  variesWithinAnyBlock = go []

  -- Two samples of one block with different values.
  Clash : List (A × B) → Set
  Clash xs = ∃₂ λ x y → x ∈ xs × y ∈ xs × proj₂ x ≡ proj₂ y × proj₁ x ≢ proj₁ y

  -- Every two samples of one block have the same value.
  Constant : List (A × B) → Set
  Constant xs = ∀ {x y} → x ∈ xs → y ∈ xs → proj₂ x ≡ proj₂ y → proj₁ x ≡ proj₁ y

  Clash⇒¬Constant : ∀ {xs} → Clash xs → ¬ Constant xs
  Clash⇒¬Constant (_ , _ , x∈ , y∈ , same , diff) c = diff (c x∈ y∈ same)

  -- The loop invariant: `seen` is the dictionary of the prefix `P` already
  -- read, and that prefix has no clash (otherwise the loop has returned).
  record Represents (seen : List (B × A)) (P : List (A × B)) : Set where
    field
      constant : Constant P
      complete : ∀ {x} → x ∈ P → firstSeen seen (proj₂ x) ≡ just (proj₁ x)
      sound    : ∀ {b w} → firstSeen seen b ≡ just w → (w , b) ∈ P

  open Represents

  start : Represents [] []
  start = record { constant = λ () ; complete = λ () ; sound = λ () }

  private
    firstSeen-hit : ∀ seen b v → firstSeen ((b , v) ∷ seen) b ≡ just v
    firstSeen-hit seen b v with b ≟ᴮ b
    ... | yes _ = refl
    ... | no b≢b = contradiction refl b≢b

    firstSeen-miss : ∀ seen {b b′} v → b′ ≢ b →
                     firstSeen ((b , v) ∷ seen) b′ ≡ firstSeen seen b′
    firstSeen-miss seen {b} {b′} v b′≢b with b′ ≟ᴮ b
    ... | yes b′≡b = contradiction b′≡b b′≢b
    ... | no _ = refl

    just-injective : ∀ {x y : A} → just x ≡ just y → x ≡ y
    just-injective refl = refl

    -- A block the dictionary has never seen has no sample in the prefix.
    unseen : ∀ {seen P b x} → Represents seen P → firstSeen seen b ≡ nothing →
             x ∈ P → proj₂ x ≢ b
    unseen rep none x∈ refl with trans (sym none) (complete rep x∈)
    ... | ()

    -- A sample read on the `nothing` branch is stored.
    extend-new : ∀ {seen P b v} → Represents seen P → firstSeen seen b ≡ nothing →
                 Represents ((b , v) ∷ seen) (P ++ [ (v , b) ])
    extend-new {seen} {P} {b} {v} rep none = record
      { constant = constant′ ; complete = complete′ ; sound = sound′ }
      where
      constant′ : Constant (P ++ [ (v , b) ])
      constant′ x∈ y∈ same with ∈-++⁻ P x∈ | ∈-++⁻ P y∈
      ... | inj₁ x∈P | inj₁ y∈P = constant rep x∈P y∈P same
      ... | inj₁ x∈P | inj₂ (here refl) = contradiction same (unseen rep none x∈P)
      ... | inj₂ (here refl) | inj₁ y∈P = contradiction (sym same) (unseen rep none y∈P)
      ... | inj₂ (here refl) | inj₂ (here refl) = refl

      complete′ : ∀ {x} → x ∈ P ++ [ (v , b) ] →
                  firstSeen ((b , v) ∷ seen) (proj₂ x) ≡ just (proj₁ x)
      complete′ x∈ with ∈-++⁻ P x∈
      ... | inj₁ x∈P = trans (firstSeen-miss seen v (unseen rep none x∈P)) (complete rep x∈P)
      ... | inj₂ (here refl) = firstSeen-hit seen b v

      sound′ : ∀ {b′ w} → firstSeen ((b , v) ∷ seen) b′ ≡ just w →
               (w , b′) ∈ P ++ [ (v , b) ]
      sound′ {b′} found with b′ ≟ᴮ b
      ... | yes refl with found
      ...   | refl = ∈-++⁺ʳ P (here refl)
      sound′ {b′} found | no _ = ∈-++⁺ˡ (sound rep found)

    -- A sample that matches its block's stored value leaves the dictionary
    -- as it was.
    extend-old : ∀ {seen P b v} → Represents seen P → firstSeen seen b ≡ just v →
                 Represents seen (P ++ [ (v , b) ])
    extend-old {seen} {P} {b} {v} rep found = record
      { constant = constant′ ; complete = complete′ ; sound = λ f → ∈-++⁺ˡ (sound rep f) }
      where
      stored : ∀ {x} → x ∈ P → proj₂ x ≡ b → proj₁ x ≡ v
      stored {x} x∈P refl = just-injective (trans (sym (complete rep x∈P)) found)

      constant′ : Constant (P ++ [ (v , b) ])
      constant′ x∈ y∈ same with ∈-++⁻ P x∈ | ∈-++⁻ P y∈
      ... | inj₁ x∈P | inj₁ y∈P = constant rep x∈P y∈P same
      ... | inj₁ x∈P | inj₂ (here refl) = stored x∈P same
      ... | inj₂ (here refl) | inj₁ y∈P = sym (stored y∈P (sym same))
      ... | inj₂ (here refl) | inj₂ (here refl) = refl

      complete′ : ∀ {x} → x ∈ P ++ [ (v , b) ] → firstSeen seen (proj₂ x) ≡ just (proj₁ x)
      complete′ x∈ with ∈-++⁻ P x∈
      ... | inj₁ x∈P = complete rep x∈P
      ... | inj₂ (here refl) = found

  -- `return true` is only reached through a clash.
  go-true : ∀ {seen P} xs → Represents seen P → go seen xs ≡ true → Clash (P ++ xs)
  go-true [] rep ()
  go-true {seen} {P} ((v , b) ∷ xs) rep h with firstSeen seen b in found
  ... | nothing =
    subst Clash (++-assoc P [ (v , b) ] xs) (go-true xs (extend-new rep found) h)
  ... | just w with w ≟ᴬ v
  ...   | yes refl =
    subst Clash (++-assoc P [ (w , b) ] xs) (go-true xs (extend-old rep found) h)
  ...   | no w≢v =
    (w , b) , (v , b) , ∈-++⁺ˡ (sound rep found) , ∈-++⁺ʳ P (here refl) , refl , w≢v

  -- `false` at the end of the loop means no two samples of a block differ.
  go-false : ∀ {seen P} xs → Represents seen P → go seen xs ≡ false → Constant (P ++ xs)
  go-false {P = P} [] rep _ = subst Constant (sym (++-identityʳ P)) (constant rep)
  go-false {seen} {P} ((v , b) ∷ xs) rep h with firstSeen seen b in found
  ... | nothing =
    subst Constant (++-assoc P [ (v , b) ] xs) (go-false xs (extend-new rep found) h)
  ... | just w with w ≟ᴬ v
  ...   | yes refl =
    subst Constant (++-assoc P [ (w , b) ] xs) (go-false xs (extend-old rep found) h)
  ...   | no _ with h
  ...     | ()

  varies-true⇒Clash : ∀ xs → variesWithinAnyBlock xs ≡ true → Clash xs
  varies-true⇒Clash xs = go-true xs start

  varies-false⇔Constant : ∀ xs → variesWithinAnyBlock xs ≡ false ⇔ Constant xs
  varies-false⇔Constant xs = mk⇔ (go-false xs start) from
    where
    from : Constant xs → variesWithinAnyBlock xs ≡ false
    from c with variesWithinAnyBlock xs in v
    ... | true = ⊥-elim (Clash⇒¬Constant (varies-true⇒Clash xs v) c)
    ... | false = refl

  ----------------------------------------------------------------------
  -- From the list Julia reads to the samples the permutations move

  -- Sample i of `xs` is `lookup xs i`: its value and its block.
  labelOf : (xs : List (A × B)) → Fin (length xs) → A
  labelOf xs i = proj₁ (lookup xs i)

  blockOf : (xs : List (A × B)) → Fin (length xs) → B
  blockOf xs i = proj₂ (lookup xs i)

  module _ (xs : List (A × B)) where
    open Labels (blockOf xs) (labelOf xs)

    Constant⇔BlockConstant : Constant xs ⇔ BlockConstant
    Constant⇔BlockConstant = mk⇔ to from
      where
      to : Constant xs → BlockConstant
      to c i j same = c (∈-lookup i) (∈-lookup j) same

      from : BlockConstant → Constant xs
      from bc x∈ y∈ same =
        let ex = lookup-index x∈
            ey = lookup-index y∈
        in trans (cong proj₁ ex)
             (trans (bc (index x∈) (index y∈)
                        (trans (sym (cong proj₂ ex)) (trans same (cong proj₂ ey))))
                    (sym (cong proj₁ ey)))

    -- Julia withholds a term's p-value (`_varies_within_any_block` is false)
    -- exactly when no within-block permutation exchanges its labels.
    withheld⇔never-exchanged :
      variesWithinAnyBlock xs ≡ false ⇔ (∀ π → WithinBlock π → Fixes π)
    withheld⇔never-exchanged = mk⇔
      (λ v → Equivalence.to untestable⇔fixed
               (Equivalence.to Constant⇔BlockConstant
                 (Equivalence.to (varies-false⇔Constant xs) v)))
      (λ fixed → Equivalence.from (varies-false⇔Constant xs)
                   (Equivalence.from Constant⇔BlockConstant
                     (Equivalence.from untestable⇔fixed fixed)))

    -- The block-constancy filter admits a term exactly when some within-block
    -- permutation exchanges its labels. Other checks can still withhold it.
    tested⇔exchanged :
      variesWithinAnyBlock xs ≡ true ⇔ Σ (Permutation′ (length xs)) λ π → WithinBlock π × ¬ Fixes π
    tested⇔exchanged = mk⇔ to from
      where
      to : variesWithinAnyBlock xs ≡ true →
           Σ (Permutation′ (length xs)) λ π → WithinBlock π × ¬ Fixes π
      to v with varies-true⇒Clash xs v
      ... | x , y , x∈ , y∈ , same , diff =
        let ex = lookup-index x∈
            ey = lookup-index y∈
        in exchanged (trans (sym (cong proj₂ ex)) (trans same (cong proj₂ ey)))
                     (λ eq → diff (trans (cong proj₁ ex) (trans eq (sym (cong proj₁ ey)))))

      from : (Σ (Permutation′ (length xs)) λ π → WithinBlock π × ¬ Fixes π) →
             variesWithinAnyBlock xs ≡ true
      from (π , wb , moved) with variesWithinAnyBlock xs in v
      ... | true = refl
      ... | false = contradiction (Equivalence.to withheld⇔never-exchanged v π wb) moved

------------------------------------------------------------------------
-- 4. The fixtures in test/unit/test_permanova.jl

module Fixtures where
  open Varies ℕₚ._≟_ ℕₚ._≟_

  -- "_varies_within_any_block": blocks i1 i1 i2 i2.
  constant-in-blocks : variesWithinAnyBlock ((1 , 1) ∷ (1 , 1) ∷ (2 , 2) ∷ (2 , 2) ∷ []) ≡ false
  constant-in-blocks = refl

  varies-in-a-block : variesWithinAnyBlock ((1 , 1) ∷ (2 , 1) ∷ (1 , 2) ∷ (1 , 2) ∷ []) ≡ true
  varies-in-a-block = refl

  -- The PERMANOVA fixture: six individuals (the blocks), each sampled in runs 1
  -- and 2; individuals 1-3 are group 1, individuals 4-6 group 2.
  individual : List ℕ
  individual = 1 ∷ 1 ∷ 2 ∷ 2 ∷ 3 ∷ 3 ∷ 4 ∷ 4 ∷ 5 ∷ 5 ∷ 6 ∷ 6 ∷ []

  pair : List ℕ → List ℕ → List (ℕ × ℕ)
  pair (v ∷ vs) (b ∷ bs) = (v , b) ∷ pair vs bs
  pair _ _ = []

  group run : List (ℕ × ℕ)
  group = pair (1 ∷ 1 ∷ 1 ∷ 1 ∷ 1 ∷ 1 ∷ 2 ∷ 2 ∷ 2 ∷ 2 ∷ 2 ∷ 2 ∷ []) individual
  run   = pair (1 ∷ 2 ∷ 1 ∷ 2 ∷ 1 ∷ 2 ∷ 1 ∷ 2 ∷ 1 ∷ 2 ∷ 1 ∷ 2 ∷ []) individual

  -- group's p-value is withheld, and rightly: no within-block permutation of
  -- the twelve samples moves a single group label.
  fixture-group-withheld :
    ∀ π → Labels.WithinBlock (blockOf group) (labelOf group) π →
          Labels.Fixes (blockOf group) (labelOf group) π
  fixture-group-withheld = Equivalence.to (withheld⇔never-exchanged group) refl

  -- run passes the block-constancy filter: some within-block permutation
  -- exchanges run labels.
  fixture-run-tested :
    Σ (Permutation′ 12) λ π → Labels.WithinBlock (blockOf run) (labelOf run) π ×
                              ¬ Labels.Fixes (blockOf run) (labelOf run) π
  fixture-run-tested = Equivalence.to (tested⇔exchanged run) refl
