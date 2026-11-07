# SPDX-FileCopyrightText: 2026 Joshua Benjamin Jewell
# SPDX-License-Identifier: AGPL-3.0-or-later

# Validation of the taxonomy source (VSEARCH or DADA2) named in a route.
using JSON3, CSV, DataFrames, OrderedCollections, DuckDB, DBInterface, Dates
using ..Categories

_validate_source(source::String) = source in ("VSEARCH", "DADA2")

function _require_study_run_source(study::String, run::String, source::String)
    err = _validate_run_request(study, run)
    isnothing(err) || return err
    _validate_source(source) || return json_error(400, "invalid_source",
                                                  "Source must be 'VSEARCH' or 'DADA2'")
    nothing
end
