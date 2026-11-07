# SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

root = joinpath(@__DIR__, "..", "..")
include(joinpath(root, "src/core/types.jl"))
include(joinpath(root, "src/core/databases.jl"))
using .Databases
Databases.ensure_databases(joinpath(root, "config/ci/databases.yml"))
