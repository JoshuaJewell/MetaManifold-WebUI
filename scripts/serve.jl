# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Starts the web server, reading JULIA_METAMANIFOLD_ROOT and JULIA_METAMANIFOLD_PORT.
#   julia --project=. scripts/serve.jl
using MetaManifold
MetaManifold.Server.main()
