{-# OPTIONS --safe --without-K #-}
-- SPDX-License-Identifier: AGPL-3.0-only
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- Must be REFUSED: reading a serial run's peak as the sum of the two peaks,
-- i.e. occupancy as cost (occupancy-types R-6).
module SerialIsCost where

open import ResourceAlgebra.Prelude
open import ResourceAlgebra.Instances.HWM
open import TypedPeak

wrong : peak (TH.grade (exits 9204 TH.∷ exits 9180 TH.∷ TH.[])) ≡ 18384
wrong = refl
