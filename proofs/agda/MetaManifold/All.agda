-- SPDX-License-Identifier: AGPL-3.0-only
-- SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
--
-- The gate entry point.  Type-checking this one module type-checks every proof
-- in the tree, which is what `just proofs` and `.github/workflows/proofs.yml`
-- do.  A module that is not imported here is not checked, so adding a module
-- without adding it here is a silent gap: `proofs/tests/axiom-audit.sh`
-- reports every `.agda` file under this directory that `All.agda` does not
-- reach, and its control proves that the report can be non-empty.

{-# OPTIONS --without-K --safe #-}

module MetaManifold.All where

-- Issue #1: Validated statistics layer — numeric core formal verification.
import MetaManifold.Prelude
import MetaManifold.Proportions
import MetaManifold.ExactCounts
import MetaManifold.PermutationTest
import MetaManifold.BenjaminiHochberg
import MetaManifold.DecimalRounding

-- Shared compositional vocabulary (also for issue #21).
import MetaManifold.Composition.Tree
import MetaManifold.Composition.Node
import MetaManifold.Composition.Sum

-- Issue #20: ILR bases from trees, SBPs and dendrograms.
import MetaManifold.ILR.SBP
import MetaManifold.ILR.Contrast
import MetaManifold.ILR.Kernel
import MetaManifold.ILR.Invariance
import MetaManifold.ILR.Orthonormal
import MetaManifold.ILR.Comb
import MetaManifold.ILR.Integer

-- Issue #7: Evidence Mode — candidate semantics, the signed finite model,
-- and the decision procedures the server routes and residual explorer run.
import MetaManifold.Evidence.Prelude
import MetaManifold.Evidence.Residual
import MetaManifold.Evidence.Echo
import MetaManifold.Evidence.Warrant
import MetaManifold.Evidence.Signed
import MetaManifold.Evidence.Decision
