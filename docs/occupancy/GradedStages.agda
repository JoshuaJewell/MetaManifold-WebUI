{-# OPTIONS --safe --without-K #-}
-- SPDX-License-Identifier: AGPL-3.0-only
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- The vsearch stage's grades, and occupancy-types' reclamation theorems at the
-- measured numbers. Fork issue #105 asks for exactly this: grade the stages
-- (alloc/hold/release) and compare the predicted peak with the measured one.
-- TypedPeak.agda computes the typed peak by the algebra's sequence law; this
-- module adds the grades the measurement itself gives, and the laws that say
-- what the difference between a peak and a cost reading is: occupancy≤cost, and
-- the strict reclamation at a cut where something is released.
--
-- Consumed, never re-declared, from the siblings scripts/fetch-siblings.sh pins
-- and scripts/check-occupancy.sh checks:
--   occupancy-types  af957280075699c659de0abd9a3b0defade77309
--   resource-algebra cf9235c532712a783a404c359a1055f40134794a
-- Everything is stated in this file's own Timeline instance (`TH`, as
-- TypedPeak.agda does): the Grade-level laws are Reclaim's as they stand, and
-- the two protocol-level ones are Reclaim's own two-step — its `cut-factor` over
-- its `reclaim-law`/`occupancy≤cost` — applied to `xs` and `ys` of that instance.
--
-- A grade is `hwm peak drift _`: what an interaction reached, and what it still
-- holds when the next one begins. The peaks below are the 2026-10-14 run's
-- (docs/occupancy-2026-10.adoc, docs/occupancy/results-2026-10-14.yml), in KiB.
module GradedStages where

open import ResourceAlgebra.Prelude
open import ResourceAlgebra.Instances.HWM
open import OccupancyTypes.Order
open import OccupancyTypes.Reclaim

import OccupancyTypes.Timeline
module TH = OccupancyTypes.Timeline HWM

-- A search that has not been waited for when the next one starts: it still holds
-- everything it reached, so the cut reclaims none of it.
still : Nat → Grade
still p = hwm p p (≤-refl p)

-- A search that has exited before the next one starts: it holds nothing, so the
-- cut reclaims all of it.
gone : Nat → Grade
gone p = hwm p 0 z≤n

-- The full pipeline runs the two searches together: ASV is still resident when
-- OTU starts, so the stage peak is their sum and nothing is reclaimed —
-- Reclaim.no-reclaim's equality, whose `refl` argument discharges the theorem's
-- own hypothesis (drift = peak) at the measured grade.
run-together : peak (still 9176 ⊗h gone 9204) ≡ 18380
run-together = no-reclaim (still 9176) refl (gone 9204)

-- A single-stage re-run runs them one after the other: ASV has exited, so the
-- cut reclaims its 9204 KiB and the stage peak is the larger search's
-- (Reclaim.reclaim-law, by `refl`).
run-in-turn : peak (gone 9204 ⊗h gone 9180) ≡ 9204
run-in-turn = refl

-- The same claim as Reclaim.reclaimed states it: with drift < peak at the cut
-- and a suffix that does anything, occupancy is strictly below cost. 9204 KiB
-- are released at the frontier, so the cost reading (18384 KiB) overstates the
-- stage.
run-in-turn<cost : peak (gone 9204 ⊗h gone 9180) < peak (gone 9204) + peak (gone 9180)
run-in-turn<cost = reclaimed (gone 9204) (gone 9180) (s≤s z≤n) (s≤s z≤n)

-- The stage's own timeline: ASV 9204 then OTU 9180, the shape of
-- docs/occupancy/vsearch-asv-then-otu.sir.
vsearch-timeline : TH.Timeline
vsearch-timeline = gone 9204 TH.∷ gone 9180 TH.∷ TH.[]

-- Prove the Reclaim.occupancy-at-cut statement for the ASV/OTU cut: the stage
-- peak is the maximum of the prefix peak and prefix drift plus suffix peak,
-- with all memory figures in KiB.
at-the-cut : peak (TH.grade vsearch-timeline) ≡
             max (peak (TH.grade (gone 9204 TH.∷ TH.[])))
                 (drift (TH.grade (gone 9204 TH.∷ TH.[])) +
                  peak (TH.grade (gone 9180 TH.∷ TH.[])))
at-the-cut =
  trans (cong peak (TH.cut-factor (gone 9204 TH.∷ TH.[]) (gone 9180 TH.∷ TH.[])))
        (reclaim-law (TH.grade (gone 9204 TH.∷ TH.[]))
                     (TH.grade (gone 9180 TH.∷ TH.[])))

-- And occupancy never exceeds cost, at every cut of a pipeline: the slack is
-- what the cut reclaims (the statement Reclaim.protocol-occupancy≤cost proves,
-- from Reclaim.occupancy≤cost along Reclaim.cut-factor).
any-cut≤cost : ∀ (xs ys : TH.Timeline)
             → peak (TH.grade (TH._++_ xs ys)) ≤ peak (TH.grade xs) + peak (TH.grade ys)
any-cut≤cost xs ys =
  subst (_≤ peak (TH.grade xs) + peak (TH.grade ys))
        (sym (cong peak (TH.cut-factor xs ys)))
        (occupancy≤cost (TH.grade xs) (TH.grade ys))

-- For the two searches as measured, the sum bound and the reclaiming bound
-- together bracket the concurrent stage's peak: 18380 KiB, no more than their
-- cost and no less than the larger search alone.
together≤cost : peak (still 9176 ⊗h gone 9204) ≤ peak (still 9176) + peak (gone 9204)
together≤cost = occupancy≤cost (still 9176) (gone 9204)
