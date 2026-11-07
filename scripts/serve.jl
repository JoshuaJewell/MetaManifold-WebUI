# SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
# SPDX-License-Identifier: AGPL-3.0-or-later

# Starts the web server, reading JULIA_METAMANIFOLD_ROOT and JULIA_METAMANIFOLD_PORT.
#   julia --project=. scripts/serve.jl
using MetaManifold
MetaManifold.Server.main()
