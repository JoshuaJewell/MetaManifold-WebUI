# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
# SPDX-License-Identifier: MPL-2.0

# The subprocess needs the production runner and logger, but no embedded R or
# web server. Include their actual sources, with no mocked protocol operations.
module RemoteProtocol
include("../../src/core/types.jl")
include("../../src/core/log.jl")
module RemoteExec
using Logging
using ..PipelineLog
include("../../src/pipeline/remote_session.jl")
end
end
