{-# OPTIONS --safe --without-K #-}
-- SPDX-License-Identifier: AGPL-3.0-only
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- Must be REFUSED: reading the concurrent stage as one that hands on nothing.
-- Running the two searches together leaves the ASV search resident when the OTU
-- search starts, so the stage's drift is that search's peak — 9176 KiB, what a
-- suffix would see — and not 0. Grading the drift as 0 would read the stage's
-- COST as its occupancy, which is the very reading the measured overlap exists
-- to refute.
module ConcurrentDriftIsZero where

open import ResourceAlgebra.Prelude
open import ResourceAlgebra.Instances.HWM
open import GradedStages

wrong : drift (still 9176 ⊗h gone 9204) ≡ 0
wrong = refl
