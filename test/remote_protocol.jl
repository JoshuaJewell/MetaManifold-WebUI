# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

# Fast, standalone protocol checks: no R installation, server or SSH account.
using Test
include("helpers/remote_protocol_module.jl")
const MetaManifold = RemoteProtocol
include("unit/remote_protocol_tests.jl")
