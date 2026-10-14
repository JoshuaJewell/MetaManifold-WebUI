{-# OPTIONS --safe --without-K #-}
-- SPDX-License-Identifier: AGPL-3.0-only
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- The vsearch stage's typed peak, computed (TESTED, by refl) with
-- occupancy-types' Timeline over resource-algebra's HWM grade, both consumed
-- by pin through scripts/fetch-siblings.sh:
--   occupancy-types  af957280075699c659de0abd9a3b0defade77309
--   resource-algebra cf9235c532712a783a404c359a1055f40134794a
-- Nothing here re-declares their definitions.
--
-- Each search is one interaction graded by its measured peak, in KiB, from the
-- 2026-10-14 run in docs/occupancy-2026-10.adoc. A search's process exits
-- holding nothing, so run alone its drift is 0. Run concurrently, the first is
-- still resident when the second starts, so that is graded with drift = peak.
-- That is the interleaving the measurement observed (both resident at once).
-- The algebra's sequence law, max (peak a) (drift a + peak b), then gives:
-- serial = the larger peak; concurrent = the sum.
module TypedPeak where

open import ResourceAlgebra.Prelude
open import ResourceAlgebra.Instances.HWM

import OccupancyTypes.Timeline
module TH = OccupancyTypes.Timeline HWM

-- A search whose process has exited before the next one starts.
exits : Nat → Grade
exits p = hwm p 0 z≤n

-- A search still resident when the next one starts.
resident : Nat → Grade
resident p = hwm p p (≤-refl p)

-- ASV 9204 then OTU 9180.
asv-then-otu : peak (TH.grade (exits 9204 TH.∷ exits 9180 TH.∷ TH.[])) ≡ 9204
asv-then-otu = refl

-- OTU 9192 then ASV 9152.
otu-then-asv : peak (TH.grade (exits 9192 TH.∷ exits 9152 TH.∷ TH.[])) ≡ 9192
otu-then-asv = refl

-- ASV 9176 and OTU 9204 together.
concurrent : peak (TH.grade (resident 9176 TH.∷ exits 9204 TH.∷ TH.[])) ≡ 18380
concurrent = refl
