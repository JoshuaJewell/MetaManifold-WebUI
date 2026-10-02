-- SPDX-License-Identifier: AGPL-3.0-only
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- The root of the proof tree.  Typechecking this file typechecks every proof
-- module; the axiom audit (proofs/tests/axiom-audit.sh) refuses any module
-- under proofs/agda/ that this file does not reach.

{-# OPTIONS --safe --without-K #-}

module MetaManifold.All where

import MetaManifold.Prelude
import MetaManifold.BenjaminiHochberg
