{-# OPTIONS --safe --without-K #-}
-- SPDX-License-Identifier: AGPL-3.0-only
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- Must be REFUSED: the claim that running the two searches concurrently leaves
-- the peak where a serial run puts it (the larger of the two).
module ConcurrentIsSerial where

open import ResourceAlgebra.Prelude
open import ResourceAlgebra.Instances.HWM
open import TypedPeak

wrong : peak (TH.grade (resident 9176 TH.∷ exits 9204 TH.∷ TH.[])) ≡ 9204
wrong = refl
