# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# Backwards compatibility shim — canonical implementation is in AnalysisConfig.jl (capital A)
# This file is kept for backwards compatibility with existing imports and CI that references analysis_config.jl
# The canonical implementation with immutable struct, validators, Nickel/DEED schemas, DANGER banner, DOI bundles is in AnalysisConfig.jl
# Fix for Julia 1.12 doc! collision: impl module is AnalysisConfigImpl, public alias is AnalysisConfig
include("AnalysisConfig.jl")

# Public alias — preserves `MetaManifold.AnalysisConfig.AnalysisConfig` API
# AnalysisConfigImpl defines struct AnalysisConfig, so no module/type name collision inside impl
const AnalysisConfig = AnalysisConfigImpl
