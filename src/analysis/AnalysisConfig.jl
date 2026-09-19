# SPDX-License-Identifier: AGPL-3.0-only
# SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
"""
    AnalysisConfig — immutable, versioned, explicit, provenance-rich analysis configuration

Implements exactly the user's answers for v1:
- Methods: NB GLM, CLR/ILR+Gaussian LM, logistic
- BH mandatory, hard-stop with DANGER banner on overrides
- Advanced Analysis section behind Evidence Mode with heavy validation/help/warnings for custom pseudocount/epsilon/zero_policy/etc.
- JSON + Nickel + DEED schemes from hyperpolymath/standards (draft 2020-12, ABNF)
- Validators that refuse meaningless inputs
- Scary DANGER banner logging for paper writers on overrides
- Full DOI-ready JSON manifest bundles with DataCite

No silent switching or auto-selection. Every field explicit, immutable.
Standards alignment:
- JSON: https://json-schema.org/draft/2020-12/schema — \$id https://hyperpolymath.github.io/MetaManifold-WebUI/schemas/analysis_config.schema.json
- Nickel: 1-formats/k9/*.ncl style, contracts ValidFormula, PseudocountContract, CorrectionContract BH mandatory, ZeroHandlingContract, MethodNormalizationCompatibility
- DEED: DEED-GRAMMAR-SPEC.adoc v0.2.0 DRAFT — :schema-version first, only () brackets, #t/#f booleans, :kebab-case keywords, filename dispatch *_chora.deed → repo-deed, SPDX header mandatory, #u5 UUID5
"""
module AnalysisConfigImpl

using Dates
using SHA
using UUIDs
using JSON3
using OrderedCollections
using Logging

# Re-use provenance from core
using ..Provenance: CapturedEnvironment, probe_metamanifold, probe_host

export AnalysisMethod, ZeroPolicy, NormalizationConfig, CorrectionConfig, AdvancedConfig, AdvancedOverrides,
       AnalysisConfig, AnalysisConfigStruct, AnalysisResult,
       validate_config, config_hash, canonical_json,
       context_help, danger_banner, is_dangerous, log_danger_banner,
       to_json, from_json, to_deed, to_nickel, from_nickel,
       create_doi_bundle, validate_nickel, validate_deed,
       AVEC_FIBRE_COLUMN, EPISTEMIC_STATUS_VALUES, DANGER_ACK_TOKEN,
       SCHEMA_VERSION, SCHEMA_VERSIONS_SUPPORTED

# --------------------------------------------------------------------------
# Constants — explicit, no silent defaults, from standards
# --------------------------------------------------------------------------

const SCHEMA_VERSION = "1.0.0"
const SCHEMA_VERSIONS_SUPPORTED = ("1.0.0",)
const AVEC_FIBRE_COLUMN = "avec_fibre"
const EPISTEMIC_STATUS_VALUES = (
    "present_in_every_admissible_world",
    "present_in_some_admissible_world",
    "absent_in_every_admissible_world",
    "unknown"
)

const DANGER_ACK_TOKEN = "I_UNDERSTAND_THE_RISK_AND_WANT_TO_OVERRIDE_BH"

@enum AnalysisMethod begin
    NB_GLM = 1        # Negative Binomial GLM (DESeq2/MASS style) — counts with overdispersion
    CLR_LM = 2        # CLR transform + Gaussian LM (Aitchison geometry)
    ILR_LM = 3        # ILR transform + Gaussian LM (balances, phylogenetic basis)
    LOGISTIC = 4      # Logistic regression (presence/absence, binary outcome)
end

const METHOD_STRINGS = Dict{String,AnalysisMethod}(
    "nb_glm" => NB_GLM, "clr_lm" => CLR_LM, "ilr_lm" => ILR_LM, "logistic" => LOGISTIC,
    "NB_GLM" => NB_GLM, "CLR_LM" => CLR_LM, "ILR_LM" => ILR_LM, "LOGISTIC" => LOGISTIC,
    "nb.glm" => NB_GLM, "clr.lm" => CLR_LM, "ilr.lm" => ILR_LM,
)

const METHOD_TO_STRING = Dict{AnalysisMethod,String}(
    NB_GLM => "nb_glm", CLR_LM => "clr_lm", ILR_LM => "ilr_lm", LOGISTIC => "logistic",
)

@enum ZeroPolicy begin
    PSEUDOCOUNT = 1
    MULTIPLICATIVE_REPLACEMENT = 2
    BAYESIAN_MULTIPLICATIVE = 3
    REFUSE = 4
end

const ZERO_POLICY_STRINGS = Dict{String,ZeroPolicy}(
    "pseudocount" => PSEUDOCOUNT,
    "multiplicative_replacement" => MULTIPLICATIVE_REPLACEMENT,
    "bayesian_multiplicative" => BAYESIAN_MULTIPLICATIVE,
    "refuse" => REFUSE,
)

const VALID_DISPERSION_METHODS = ("parametric", "local", "mean", "pooled", "glmGamPoi")
const VALID_ZERO_HANDLING = ("pseudocount", "multiplicative_replacement", "bayesian_multiplicative", "refuse")
const VALID_ILR_BASIS = ("default", "phylogenetic", "sequential_binary_partition", "balance_dendrogram")
const VALID_NORMALIZATION_FOR_METHOD = Dict{AnalysisMethod, Vector{String}}(
    NB_GLM => ["none", "rarefy", "size_factors", "relative", "TSS", "CSS", "RSS"], # TSS/CSS/RSS deferred but allowed as alias for relative
    CLR_LM => ["clr"],
    ILR_LM => ["ilr"],
    LOGISTIC => ["none", "relative", "rarefy", "presence_absence", "TSS"],
)

# --------------------------------------------------------------------------
# Sub-configs — Advanced Analysis section with heavy validation/help/warnings
# --------------------------------------------------------------------------

"""
    NormalizationConfig(; method, ...)

Create a normalisation or compositional-transform configuration. The method name
is normalised to lower case and the zero-policy name is parsed case-insensitively.
Incompatible options, non-positive CLR/ILR pseudocounts, and invalid epsilon or
replacement-delta values raise `ArgumentError`; unusual accepted values may emit
warnings.
"""
struct NormalizationConfig
    method::String
    pseudocount::Float64
    epsilon::Float64
    zero_policy::ZeroPolicy
    ilr_basis::Union{String,Nothing}
    multiplicative_replacement_delta::Union{Float64,Nothing}
    tss_css_rss_note::Union{String,Nothing}  # deferred feature note

    function NormalizationConfig(;
        method::String,
        pseudocount::Float64=0.5,
        epsilon::Float64=1e-6,
        zero_policy::String="pseudocount",
        ilr_basis::Union{String,Nothing}=nothing,
        multiplicative_replacement_delta::Union{Float64,Nothing}=nothing,
        tss_css_rss_note::Union{String,Nothing}=nothing
    )
        method_clean = lowercase(strip(method))
        isempty(method_clean) && throw(ArgumentError("normalization.method must be non-empty (e.g. 'clr', 'size_factors', 'TSS') — see context_help('normalization.method')"))

        # Refuse meaningless: ; backtick dollar injection
        if occursin(r"[;`dollar]", method_clean)
            throw(ArgumentError("normalization.method contains forbidden ; ` \$ (injection prevention) — got '\$method_clean'"))
        end

        # Zero policy parsing
        zp_str = lowercase(strip(zero_policy))
        zp = get(ZERO_POLICY_STRINGS, zp_str, nothing)
        isnothing(zp) && throw(ArgumentError("zero_policy must be one of $(join(keys(ZERO_POLICY_STRINGS), ", ")) — got '\$zero_policy'. See context_help('advanced.zero_policy')"))

        # Heavy validation per method
        if method_clean in ("clr", "ilr")
            pseudocount <= 0 && throw(ArgumentError("For compositional methods (CLR/ILR), pseudocount must be >0 (got \$pseudocount). Zero replacement mandatory because log(0) undefined. See context_help('normalization.pseudocount')"))
            pseudocount >= 1 && @warn "pseudocount >=1 unusual for CLR/ILR (got \$pseudocount); typical 0.5 or 0.65. May distort low-abundance features." pseudocount method_clean
            if zp == REFUSE
                throw(ArgumentError("zero_policy='refuse' is mathematically invalid for CLR/ILR (log(0) undefined). Refusing even with DANGER token. Use pseudocount or multiplicative_replacement."))
            end
            if method_clean == "ilr"
                if !isnothing(ilr_basis) && !(ilr_basis in VALID_ILR_BASIS)
                    throw(ArgumentError("ilr_basis must be one of $(join(VALID_ILR_BASIS, ", ")) (got '\$ilr_basis')"))
                end
            else # clr
                if !isnothing(ilr_basis)
                    throw(ArgumentError("ilr_basis is meaningless for CLR (only for ILR). Refusing. See context_help('normalization.ilr_basis')"))
                end
            end
        else
            if !isnothing(ilr_basis)
                throw(ArgumentError("ilr_basis only meaningful for ILR method, not for '\$method_clean'"))
            end
            # For NB_GLM, pseudocount is allowed but warned if used with size_factors
            if method_clean == "size_factors" && pseudocount != 0.5
                @warn "pseudocount is ignored for size_factors (NB_GLM). Using size_factors from DESeq2, not pseudocount. Got pseudocount=\$pseudocount — will be ignored unless you switch to clr/ilr."
            end
        end

        # Epsilon validation — advanced
        if !(0 < epsilon < 1)
            throw(ArgumentError("epsilon must be in (0,1) for numerical stability, got \$epsilon. Typical 1e-6. See context_help('advanced.epsilon')"))
        end
        if epsilon > 1e-3
            @warn "epsilon >1e-3 is large and may affect zero handling and log transforms" epsilon
        end

        # Multiplicative replacement delta
        if !isnothing(multiplicative_replacement_delta)
            delta = multiplicative_replacement_delta
            (delta <= 0 || delta >= 1) && throw(ArgumentError("multiplicative_replacement_delta must be in (0,1), got \$delta — see context_help('advanced.zero_policy')"))
        end

        # TSS/CSS/RSS note — deferred feature
        if method_clean in ("tss", "css", "rss")
            if isnothing(tss_css_rss_note)
                @warn "TSS/CSS/RSS are deferred features (see GitHub issues). Currently aliased to relative/TSS. For exact TSS/CSS/RSS offsets, see deferred issue with value/difficulty/risk."
            end
        end

        new(method_clean, pseudocount, epsilon, zp, ilr_basis, multiplicative_replacement_delta, tss_css_rss_note)
    end
end

"""
    CorrectionConfig(; method="BH", alpha=0.05, ...)

Create a multiple-testing correction configuration. Non-BH methods are accepted
only when `allow_no_correction` is true and `acknowledgment_token` matches
`DANGER_ACK_TOKEN`. Empty or unsafe method strings, invalid tokens, and alpha
values outside `(0, 1)` raise `ArgumentError`.
"""
struct CorrectionConfig
    method::String
    alpha::Float64
    allow_no_correction::Bool
    acknowledgment_token::Union{String,Nothing}

    function CorrectionConfig(;
        method::String="BH",
        alpha::Float64=0.05,
        allow_no_correction::Bool=false,
        acknowledgment_token::Union{String,Nothing}=nothing
    )
        method_clean = strip(method)
        isempty(method_clean) && throw(ArgumentError("correction.method must be non-empty — see context_help('correction.method')"))

        if occursin(r"[;`dollar]", method_clean)
            throw(ArgumentError("correction.method contains forbidden ; ` \$ — got '\$method_clean'"))
        end

        # BH mandatory check
        is_bh = uppercase(method_clean) in ("BH", "FDR", "BENJAMINI-HOCHBERG", "BENJAMINI_HOCHBERG")
        if !is_bh
            if !allow_no_correction
                throw(ArgumentError("p-value correction method must be BH (Benjamini-Hochberg) in v1. Got '\$method_clean'. Microbiome data tests thousands of taxa — uncorrected p-values give ~5% false positives under null. If you truly want to override, set allow_no_correction=true and acknowledgment_token='\$DANGER_ACK_TOKEN'. This will trigger DANGER banner and be recorded in provenance. See context_help('correction.method')"))
            end
        end

        if allow_no_correction
            if isnothing(acknowledgment_token) || acknowledgment_token != DANGER_ACK_TOKEN
                throw(ArgumentError("DANGER: Attempting to disable BH correction — scientifically dangerous for high-dimensional data, will inflate false discoveries. To proceed, set acknowledgment_token to exactly '\$DANGER_ACK_TOKEN'. This action will be logged, bannered, and included in DOI bundle provenance. See context_help('correction.method') and danger_banner()."))
            end
        end

        (alpha <= 0 || alpha >= 1) && throw(ArgumentError("alpha must be in (0,1), got \$alpha. Typical 0.05. See context_help('correction.alpha')"))

        canonical_method = allow_no_correction ? method_clean : "BH"

        new(canonical_method, alpha, allow_no_correction, acknowledgment_token)
    end
end

"""
    AdvancedConfig(; ...)

Create the advanced analysis settings. Invalid categorical choices and values
outside their supported ranges raise `ArgumentError`; risky but accepted values
emit warnings. Either zero-handling option set to `"refuse"` also requires the
danger acknowledgement token.
"""
struct AdvancedConfig
    dispersion_method::String
    zero_handling::String
    zero_policy::ZeroPolicy
    pseudocount::Float64
    epsilon::Float64
    min_prevalence::Float64
    min_abundance::Float64
    max_features::Union{Int,Nothing}
    min_samples_per_group::Int
    robust::Bool
    acknowledgment_token::Union{String,Nothing}

    function AdvancedConfig(;
        dispersion_method::String="parametric",
        zero_handling::String="pseudocount",
        zero_policy::String="pseudocount",
        pseudocount::Float64=0.5,
        epsilon::Float64=1e-6,
        min_prevalence::Float64=0.1,
        min_abundance::Float64=0.0,
        max_features::Union{Int,Nothing}=nothing,
        min_samples_per_group::Int=3,
        robust::Bool=false,
        acknowledgment_token::Union{String,Nothing}=nothing
    )
        dispersion_method_clean = lowercase(strip(dispersion_method))
        zero_handling_clean = lowercase(strip(zero_handling))
        zero_policy_clean = lowercase(strip(zero_policy))

        # Dispersion method validation
        if !(dispersion_method_clean in VALID_DISPERSION_METHODS)
            throw(ArgumentError("dispersion_method must be one of $(join(VALID_DISPERSION_METHODS, ", ")) — got '\$dispersion_method_clean'. See context_help('advanced.dispersion_method')"))
        end

        # Zero handling validation
        if !(zero_handling_clean in VALID_ZERO_HANDLING)
            throw(ArgumentError("zero_handling must be one of $(join(VALID_ZERO_HANDLING, ", ")) — got '\$zero_handling_clean'. See context_help('advanced.zero_handling')"))
        end

        # Zero policy enum
        zp = get(ZERO_POLICY_STRINGS, zero_policy_clean, nothing)
        isnothing(zp) && throw(ArgumentError("zero_policy must be one of $(join(keys(ZERO_POLICY_STRINGS), ", ")) — got '\$zero_policy_clean'"))

        # Pseudocount heavy validation
        if pseudocount <= 0
            throw(ArgumentError("advanced.pseudocount must be >0 (got \$pseudocount) because log(0) undefined. Typical 0.5. See context_help('normalization.pseudocount')"))
        end
        if pseudocount >= 1
            @warn "advanced.pseudocount >=1 unusual (got \$pseudocount); typical 0.5 or 0.65. May distort low-abundance features." pseudocount
        end
        if pseudocount < 0.1
            @warn "advanced.pseudocount <0.1 very small (got \$pseudocount); will create extreme log-ratios for zeros. Consider 0.5." pseudocount
        end

        # Epsilon heavy validation
        if !(0 < epsilon < 1)
            throw(ArgumentError("advanced.epsilon must be in (0,1) for numerical stability, got \$epsilon. Typical 1e-6. See context_help('advanced.epsilon')"))
        end
        if epsilon > 1e-3
            @warn "advanced.epsilon >1e-3 large (got \$epsilon) may affect zero handling and log transforms" epsilon
        end
        if epsilon < 1e-12
            @warn "advanced.epsilon <1e-12 extremely small (got \$epsilon) may cause underflow" epsilon
        end

        # Prevalence validation
        if !(0 <= min_prevalence <= 1)
            throw(ArgumentError("min_prevalence must be in [0,1], got \$min_prevalence. 0.1 = present in >=10% samples. See context_help('advanced.min_prevalence')"))
        end

        # Abundance validation
        if min_abundance < 0
            throw(ArgumentError("min_abundance must be >=0, got \$min_abundance"))
        end

        # Max features validation
        if !isnothing(max_features)
            if max_features <= 0
                throw(ArgumentError("max_features must be >0 or nothing, got \$max_features — meaningless to test 0 features"))
            end
            if max_features > 100000
                throw(ArgumentError("max_features >100000 (got \$max_features) is excessive and will cause memory issues. Refusing."))
            end
            if max_features < 10
                @warn "max_features <10 very small (got \$max_features) — will test only \$max_features features, may miss biology" max_features
            end
        end

        # Min samples per group validation
        if min_samples_per_group < 2
            throw(ArgumentError("min_samples_per_group must be >=2 (got \$min_samples_per_group) — need at least 2 samples per group for statistical test, 3 recommended"))
        end
        if min_samples_per_group < 3
            @warn "min_samples_per_group <3 (got \$min_samples_per_group) — statistical power very low, results may be unreliable" min_samples_per_group
        end

        # Zero handling refuse requires token — DANGER
        if zero_handling_clean == "refuse" || zp == REFUSE
            if isnothing(acknowledgment_token) || acknowledgment_token != DANGER_ACK_TOKEN
                throw(ArgumentError("DANGER: zero_handling='refuse' will cause log(0) for CLR/ILR and biased handling for NB_GLM. Requires acknowledgment_token='\$DANGER_ACK_TOKEN'. Even then, CLR/ILR + refuse is mathematically invalid and will be refused at runtime. See context_help('advanced.zero_handling')"))
            end
        end

        new(dispersion_method_clean, zero_handling_clean, zp, pseudocount, epsilon, min_prevalence, min_abundance, max_features, min_samples_per_group, robust, acknowledgment_token)
    end
end

# --------------------------------------------------------------------------
# Main immutable AnalysisConfig struct — exactly user's answers
# --------------------------------------------------------------------------

"""
    AnalysisConfig(; method, formula, metadata_columns, ...)

Create an immutable analysis configuration after validating its schema version,
UUID, formula, metadata columns, outcome, and method/normalisation pairing. The
constructor enriches a copy of `provenance`, computes a content hash when none is
supplied, and derives the stored danger flag from correction and advanced
settings; an explicitly supplied `dangerous=true` also marks the configuration.
Invalid combinations raise `ArgumentError`, while accepted but questionable
settings may emit warnings.
"""
struct AnalysisConfig
    schema_version::String
    id::String
    created_at::DateTime
    created_by::String
    method::AnalysisMethod
    formula::String
    outcome_column::Union{String,Nothing}
    metadata_columns::Vector{String}
    normalization::NormalizationConfig
    correction::CorrectionConfig
    advanced::AdvancedConfig
    provenance::OrderedDict{String,Any}
    hash::String
    dangerous::Bool

    function AnalysisConfig(;
        schema_version::String=SCHEMA_VERSION,
        id::String=string(uuid4()),
        created_at::DateTime=now(Dates.UTC),
        created_by::String="anonymous",
        method::String,
        formula::String,
        outcome_column::Union{String,Nothing}=nothing,
        metadata_columns::Vector{String},
        normalization::NormalizationConfig=NormalizationConfig(method="size_factors"),
        correction::CorrectionConfig=CorrectionConfig(),
        advanced::AdvancedConfig=AdvancedConfig(),
        provenance::OrderedDict{String,Any}=OrderedDict{String,Any}(),
        hash::Union{String,Nothing}=nothing,
        dangerous::Union{Bool,Nothing}=nothing
    )
        # Schema version check
        if !(schema_version in SCHEMA_VERSIONS_SUPPORTED)
            throw(ArgumentError("schema_version must be one of $(join(SCHEMA_VERSIONS_SUPPORTED, ", ")) — got '\$schema_version'. See DEED spec :schema-version first"))
        end

        # ID validation UUID
        try
            UUID(id)
        catch
            throw(ArgumentError("id must be valid UUID4, got '\$id'"))
        end

        # Method parsing
        method_clean = strip(method)
        isempty(method_clean) && throw(ArgumentError("method must be non-empty — one of nb_glm, clr_lm, ilr_lm, logistic. See context_help('method')"))
        occursin(r"[;`dollar]", method_clean) && throw(ArgumentError("method contains forbidden ; ` \$ — got '\$method_clean'"))
        method_enum = get(METHOD_STRINGS, method_clean, get(METHOD_STRINGS, lowercase(method_clean), nothing))
        isnothing(method_enum) && throw(ArgumentError("method must be one of $(join(keys(METHOD_STRINGS), ", ")) — got '\$method_clean'. No auto-selection. See context_help('method')"))

        # Formula validation — heavy, refuses meaningless
        formula_clean = strip(formula)
        isempty(formula_clean) && throw(ArgumentError("formula must be non-empty, e.g. '~ group' or 'disease ~ group + batch'. See context_help('formula')"))
        length(formula_clean) < 2 && throw(ArgumentError("formula too short, must reference at least one metadata column — got '\$formula_clean'"))
        !occursin("~", formula_clean) && throw(ArgumentError("formula must contain '~' (R-style), e.g. '~ group' — got '\$formula_clean'. See context_help('formula')"))
        if occursin(r"[;`dollar]", formula_clean)
            throw(ArgumentError("formula contains forbidden ; ` \$ (injection prevention) — got '\$formula_clean'. See context_help('formula')"))
        end
        # Refuse formulas that are just "~" or "~  "
        if strip(replace(formula_clean, "~" => "")) == ""
            throw(ArgumentError("formula must reference at least one metadata column after '~' — got '\$formula_clean'. Refusing meaningless input."))
        end

        # Metadata columns validation
        if isempty(metadata_columns)
            throw(ArgumentError("metadata_columns must be non-empty, at least 1 explicit column, no auto-selection. See context_help('metadata_columns')"))
        end
        if length(unique(metadata_columns)) != length(metadata_columns)
            throw(ArgumentError("metadata_columns must be unique, got duplicates in \$metadata_columns"))
        end
        for col in metadata_columns
            isempty(strip(col)) && throw(ArgumentError("metadata_columns contains empty string — refusing"))
            if occursin(r"[;`dollar]", col)
                throw(ArgumentError("metadata_columns contains forbidden ; ` \$ in '\$col'"))
            end
            if !occursin(r"^[a-zA-Z0-9_\.\-]+\$", col)
                throw(ArgumentError("metadata_columns must match pattern ^[a-zA-Z0-9_.\\-]+\$ — got '\$col'. See JSON schema."))
            end
        end

        # Outcome column for logistic
        if method_enum == LOGISTIC
            if isnothing(outcome_column) || isempty(strip(outcome_column))
                throw(ArgumentError("outcome_column is required for logistic regression (binary outcome). Got nothing. Formula should be 'disease ~ group' and outcome_column='disease'. See context_help('outcome_column')"))
            end
            # Outcome must be in metadata_columns or formula left side
            oc_clean = strip(outcome_column)
            # Allow outcome_column to be in metadata_columns OR left side of formula
            # For simplicity, require it in metadata_columns for now, but warn if not
            if !(oc_clean in metadata_columns)
                @warn "outcome_column '\$oc_clean' not in metadata_columns \$metadata_columns — may be left side of formula, but should be listed in metadata_columns for explicitness" outcome_column metadata_columns
            end
        else
            if !isnothing(outcome_column) && !isempty(strip(outcome_column))
                @warn "outcome_column is only used for logistic, but method is $(METHOD_TO_STRING[method_enum]) — outcome_column will be ignored unless formula uses it" outcome_column method_enum
            end
        end

        # Normalization compatibility with method
        norm_method = normalization.method
        allowed_norms = get(VALID_NORMALIZATION_FOR_METHOD, method_enum, String[])
        if !(norm_method in allowed_norms)
            throw(ArgumentError("normalization.method '\$norm_method' incompatible with method '$(METHOD_TO_STRING[method_enum])'. Allowed for $(METHOD_TO_STRING[method_enum]): $(join(allowed_norms, ", ")). See context_help('normalization.method') and MethodNormalizationCompatibility contract in Nickel. Refusing meaningless combination."))
        end

        # Dangerous computed
        is_dang = false
        if correction.allow_no_correction
            is_dang = true
        end
        if advanced.zero_handling == "refuse" || advanced.zero_policy == REFUSE
            is_dang = true
        end
        if advanced.min_samples_per_group < 3
            is_dang = true
        end
        if !isnothing(dangerous)
            is_dang = dangerous || is_dang
        end

        # Provenance enriched
        prov = OrderedDict{String,Any}(provenance)
        if !haskey(prov, "schema_version")
            prov["schema_version"] = schema_version
        end
        if !haskey(prov, "created_at")
            prov["created_at"] = string(created_at)
        end
        if !haskey(prov, "created_by")
            prov["created_by"] = created_by
        end
        if !haskey(prov, "method")
            prov["method"] = METHOD_TO_STRING[method_enum]
        end
        if !haskey(prov, "dangerous")
            prov["dangerous"] = is_dang
        end

        # Hash computation content-addressed (immutable)
        hash_input = if isnothing(hash)
            # Canonical JSON of config without hash and provenance hash to avoid circular
            canonical = OrderedDict{String,Any}(
                "schema_version" => schema_version,
                "id" => id,
                "created_at" => string(created_at),
                "created_by" => created_by,
                "method" => METHOD_TO_STRING[method_enum],
                "formula" => formula_clean,
                "outcome_column" => outcome_column,
                "metadata_columns" => metadata_columns,
                "normalization" => OrderedDict(
                    "method" => normalization.method,
                    "pseudocount" => normalization.pseudocount,
                    "epsilon" => normalization.epsilon,
                    "zero_policy" => string(normalization.zero_policy),
                    "ilr_basis" => normalization.ilr_basis,
                    "multiplicative_replacement_delta" => normalization.multiplicative_replacement_delta
                ),
                "correction" => OrderedDict(
                    "method" => correction.method,
                    "alpha" => correction.alpha,
                    "allow_no_correction" => correction.allow_no_correction
                ),
                "advanced" => OrderedDict(
                    "dispersion_method" => advanced.dispersion_method,
                    "zero_handling" => advanced.zero_handling,
                    "zero_policy" => string(advanced.zero_policy),
                    "pseudocount" => advanced.pseudocount,
                    "epsilon" => advanced.epsilon,
                    "min_prevalence" => advanced.min_prevalence,
                    "min_abundance" => advanced.min_abundance,
                    "max_features" => advanced.max_features,
                    "min_samples_per_group" => advanced.min_samples_per_group,
                    "robust" => advanced.robust
                )
            )
            bytes2hex(sha256(JSON3.write(canonical)))
        else
            hash
        end

        new(schema_version, id, created_at, created_by, method_enum, formula_clean, outcome_column, metadata_columns, normalization, correction, advanced, prov, hash_input, is_dang)
    end
end

# Backwards compatibility aliases — lowercase file and old tests use these names
const AnalysisConfigStruct = AnalysisConfig
const AdvancedOverrides = AdvancedConfig

"""
    AnalysisResult(; config_id, config_hash, method, ...)

Create an immutable result linked to an analysis configuration. The constructor
adds the configuration identifiers and method to a copy of `provenance` and
computes a content hash unless one is supplied. Invalid UUIDs or an empty
configuration hash raise `ArgumentError`.
"""
struct AnalysisResult
    id::String
    config_id::String
    config_hash::String
    created_at::DateTime
    method::AnalysisMethod
    results::OrderedDict{String,Any}
    provenance::OrderedDict{String,Any}
    hash::String

    function AnalysisResult(;
        id::String=string(uuid4()),
        config_id::String,
        config_hash::String,
        created_at::DateTime=now(Dates.UTC),
        method::AnalysisMethod,
        results::OrderedDict{String,Any}=OrderedDict{String,Any}(),
        provenance::OrderedDict{String,Any}=OrderedDict{String,Any}(),
        hash::Union{String,Nothing}=nothing
    )
        try UUID(id) catch; throw(ArgumentError("id must be valid UUID4")) end
        try UUID(config_id) catch; throw(ArgumentError("config_id must be valid UUID4")) end
        isempty(config_hash) && throw(ArgumentError("config_hash must be non-empty"))

        prov = OrderedDict{String,Any}(provenance)
        prov["config_id"] = config_id
        prov["config_hash"] = config_hash
        prov["method"] = METHOD_TO_STRING[method]

        hash_computed = if isnothing(hash)
            canonical = OrderedDict(
                "id" => id,
                "config_id" => config_id,
                "config_hash" => config_hash,
                "created_at" => string(created_at),
                "method" => METHOD_TO_STRING[method],
                "results" => results
            )
            bytes2hex(sha256(JSON3.write(canonical)))
        else
            hash
        end

        new(id, config_id, config_hash, created_at, method, results, prov, hash_computed)
    end
end

# --------------------------------------------------------------------------
# Validators — refuse meaningless inputs
# --------------------------------------------------------------------------

"""
    validate_config(config, available_columns; strict=true)

Return validation messages for missing metadata columns and formula references.
With `strict=true`, also recheck the method/normalisation pairing. A dangerous
configuration emits a warning but does not add a validation message.
"""
function validate_config(config::AnalysisConfig, available_columns::Vector{String}; strict::Bool=true)
    errors = String[]

    # Check metadata_columns exist in available
    for col in config.metadata_columns
        if !(col in available_columns)
            push!(errors, "metadata_columns '\$col' not found in available columns $(available_columns) — see context_help('metadata_columns')")
        end
    end

    # Check formula references only metadata_columns
    # Simple parsing: extract tokens after ~ and split by + * : etc.
    formula_body = config.formula
    # Remove left side if logistic
    if occursin("~", formula_body)
        parts = split(formula_body, "~")
        if length(parts) == 2
            left = strip(parts[1])
            right = strip(parts[2])
            # Left side for logistic should be outcome_column if present
            if config.method == LOGISTIC && !isempty(left)
                if !isnothing(config.outcome_column) && left != config.outcome_column
                    push!(errors, "For logistic, left side of formula '\$left' should match outcome_column '$(config.outcome_column)' — see context_help('formula')")
                end
            end
            # Right side tokens
            # Split by + * : | / etc.
            tokens = split(right, r"[\+\*\:\|\(\) ]+")
            for tok in tokens
                tok_clean = strip(tok)
                isempty(tok_clean) && continue
                # Skip R functions like I, 1, etc.
                if tok_clean in ("I", "1", "0", "")
                    continue
                end
                # Skip numeric
                if occursin(r"^\d+\$", tok_clean)
                    continue
                end
                if !(tok_clean in config.metadata_columns)
                    push!(errors, "Formula references '\$tok_clean' which is not in metadata_columns $(config.metadata_columns) — refusing. See context_help('formula')")
                end
            end
        end
    end

    # Check normalization compatibility already done in constructor, but re-check for strict
    if strict
        norm_method = config.normalization.method
        allowed = get(VALID_NORMALIZATION_FOR_METHOD, config.method, String[])
        if !(norm_method in allowed)
            push!(errors, "Incompatible normalization.method '\$norm_method' for method '$(METHOD_TO_STRING[config.method])' — allowed $(join(allowed, ", "))")
        end
    end

    # Check dangerous
    if config.dangerous
        @warn "Config is marked dangerous — will trigger DANGER banner" config_id=config.id dangerous=config.dangerous
    end

    return errors
end

"""Return the content hash stored on `config`."""
function config_hash(config::AnalysisConfig)
    return config.hash
end

"""
    canonical_json(config)

Serialise selected fields from `config` in a fixed order. Unlike
[`to_json`](@ref), this representation omits provenance and several advanced and
normalisation fields.
"""
function canonical_json(config::AnalysisConfig)
    return JSON3.write(OrderedDict(
        "schema_version" => config.schema_version,
        "id" => config.id,
        "created_at" => string(config.created_at),
        "created_by" => config.created_by,
        "method" => METHOD_TO_STRING[config.method],
        "formula" => config.formula,
        "outcome_column" => config.outcome_column,
        "metadata_columns" => config.metadata_columns,
        "normalization" => OrderedDict(
            "method" => config.normalization.method,
            "pseudocount" => config.normalization.pseudocount,
            "epsilon" => config.normalization.epsilon,
            "zero_policy" => string(config.normalization.zero_policy),
            "ilr_basis" => config.normalization.ilr_basis
        ),
        "correction" => OrderedDict(
            "method" => config.correction.method,
            "alpha" => config.correction.alpha,
            "allow_no_correction" => config.correction.allow_no_correction
        ),
        "advanced" => OrderedDict(
            "dispersion_method" => config.advanced.dispersion_method,
            "zero_handling" => config.advanced.zero_handling,
            "pseudocount" => config.advanced.pseudocount,
            "epsilon" => config.advanced.epsilon,
            "min_prevalence" => config.advanced.min_prevalence
        ),
        "hash" => config.hash,
        "dangerous" => config.dangerous
    ))
end

# --------------------------------------------------------------------------
# Context-sensitive help — for UI
# --------------------------------------------------------------------------

"""
    context_help(field_path)

Return the UI help text for a recognised analysis field path. Unknown paths
return a fallback message that includes the requested path.
"""
function context_help(field_path::String)
    help_db = Dict{String,String}(
        "method" => """
        Analysis Method (required, explicit, no auto-selection) — v1: NB GLM, CLR/ILR+Gaussian, logistic

        - nb_glm: Negative Binomial GLM for raw counts with overdispersion. Uses DESeq2-style size_factors. Best for counts, handles library size via size_factors, dispersion via parametric/local/mean/pooled/glmGamPoi. See Love et al. 2014.
        - clr_lm: Centered Log-Ratio + Gaussian LM (compositional, Aitchison geometry). Requires pseudocount >0 because log(0) undefined. Handles compositional data (relative abundances). See Gloor et al. 2017.
        - ilr_lm: Isometric Log-Ratio + Gaussian LM (balances, phylogenetic basis possible). Requires pseudocount >0 and ilr_basis. Basis options: default, phylogenetic, sequential_binary_partition, balance_dendrogram. See Egozcue et al. 2003.
        - logistic: Logistic regression for binary outcome (presence/absence). Requires outcome_column. Normalization presence_absence or relative.

        No silent switching. Every analysis explicit. See JSON schema and Nickel contract MethodNormalizationCompatibility.
        """,
        "formula" => """
        R-style formula, e.g. '~ group' or 'disease ~ group + batch'

        - Left of ~ is outcome (required for logistic, must match outcome_column)
        - Right of ~ lists metadata columns, e.g. group + batch
        - Must reference only columns in metadata_columns — no auto-selection
        - Forbidden: ; ` \$ (injection prevention, see ValidFormula contract in Nickel)
        - Must contain ~ and at least one column after ~
        - Example: "~ group" tests effect of group, "~ group + batch" controls for batch
        - For logistic: "disease ~ group" with outcome_column="disease"

        Refuses meaningless inputs: empty, "~", "group" without ~, forbidden chars.
        See JSON schema pattern ^[^;`\$]+\$ and Nickel ValidFormula.
        """,
        "outcome_column" => """
        Binary outcome column for logistic regression (required for logistic, ignored otherwise)

        - For LOGISTIC: must be non-empty and should be in metadata_columns for explicitness
        - Should match left side of formula, e.g. formula "disease ~ group" → outcome_column "disease"
        - Must be binary (0/1 or true/false) in actual data — validated at runtime
        - For NB_GLM/CLR_LM/ILR_LM: ignored unless formula uses it, but warning issued

        See context_help('formula') and JSON schema.
        """,
        "metadata_columns" => """
        Explicit list of metadata columns used — no auto-selection, must exist in study metadata

        - Min 1, unique, pattern ^[a-zA-Z0-9_.\\-]+\$ (alphanumeric + _ . -)
        - Must exist in available_columns at validation time (see validate_config)
        - No silent switching — every column explicit
        - Example: ["group", "batch", "age"] — then formula can use group + batch, but not age unless listed

        Refuses: empty list, duplicates, empty strings, forbidden ; ` \$.
        See JSON schema and Nickel contract.
        """,
        "normalization.method" => """
        Normalization / Transform (method-dependent) — must be compatible with method

        - For NB_GLM: none, size_factors (DESeq2 default, preferred), relative, rarefy (discouraged, use with caution), TSS (alias for relative, deferred exact TSS), CSS (deferred), RSS (deferred)
        - For CLR_LM: must be clr — Centered Log-Ratio, requires pseudocount >0
        - For ILR_LM: must be ilr — Isometric Log-Ratio, requires pseudocount >0 and ilr_basis
        - For LOGISTIC: presence_absence, none, relative, rarefy, TSS

        TSS/CSS/RSS offsets are deferred features (see GitHub issues) — currently aliased to relative. For exact TSS/CSS/RSS offsets, see deferred issue with value/difficulty/risk.

        Refuses meaningless: NB_GLM + clr/ilr (counts vs compositional), CLR_LM + none, etc. See Nickel MethodNormalizationCompatibility contract.
        See JSON schema enum and DEED (normalization :method).
        """,
        "normalization.pseudocount" => """
        Pseudocount for zero replacement (CLR/ILR mandatory, NB_GLM optional but warned)

        - Must be >0 because log(0) undefined — refuses 0 or negative
        - Typical: 0.5 (common), 0.65 (Martín-Fernández et al.), 1.0 (conservative but distorts)
        - Too small (<0.1) creates extreme log-ratios for zeros, too large (>=1) distorts low-abundance features — warnings issued
        - For NB_GLM size_factors: ignored, warning issued if not default 0.5
        - Advanced: custom pseudocount behind Advanced Analysis expander, hidden unless Evidence Mode, heavy validation, warnings

        See JSON schema exclusiveMinimum 0 and Nickel PseudocountContract.
        See context_help('advanced.pseudocount') for advanced warnings.
        """,
        "normalization.epsilon" => """
        Epsilon for numerical stability (advanced, behind Advanced Analysis)

        - Must be in (0,1), typical 1e-6
        - Used for zero handling and log transforms numerical stability
        - Too large >1e-3 may affect zero handling and log transforms — warning
        - Too small <1e-12 may cause underflow — warning
        - Heavy validation, refusal if not in (0,1)

        See context_help('advanced.epsilon') and Nickel contract.
        """,
        "normalization.zero_policy" => """
        Zero handling policy (advanced)

        - pseudocount: add pseudocount to zeros (default, safe)
        - multiplicative_replacement: replace zeros via multiplicative replacement (Martín-Fernández et al.), requires multiplicative_replacement_delta in (0,1)
        - bayesian_multiplicative: Bayesian multiplicative replacement
        - refuse: refuse to handle zeros — DANGEROUS, requires DANGER token, mathematically invalid for CLR/ILR (log(0) undefined), will be refused at runtime even with token

        See context_help('advanced.zero_policy') and Nickel ZeroHandlingContract.
        """,
        "normalization.ilr_basis" => """
        ILR basis (only for ILR, meaningless otherwise)

        - default: default ILR basis from compositions package
        - phylogenetic: phylogenetic tree based balances (requires tree, deferred)
        - sequential_binary_partition: SBP from user-provided partition (deferred)
        - balance_dendrogram: balance dendrogram

        Refuses meaningless use for non-ILR methods. See JSON schema enum and Nickel.
        """,
        "correction.method" => """
        Multiple testing correction — BH mandatory in v1, hard-stop DANGER banner on overrides

        Microbiome data tests thousands of taxa. Uncorrected p-values give ~5% false positives under null (e.g., 1500 taxa → 75 false positives). BH (Benjamini-Hochberg FDR) controls false discovery rate.

        - BH mandatory in v1 — any override triggers DANGER banner and requires acknowledgment token
        - Allowed: BH, FDR, Benjamini-Hochberg (all normalized to BH)
        - Override: set allow_no_correction=true and acknowledgment_token='I_UNDERSTAND_THE_RISK_AND_WANT_TO_OVERRIDE_BH' — will be logged, bannered, included in DOI bundle provenance, scary DANGER banner for paper writers

        See JSON schema and Nickel CorrectionContract, and danger_banner().

        Scientific value: Prevents p-hacking and false discoveries. See Benjamini & Hochberg 1995.
        """,
        "correction.alpha" => """
        FDR alpha (significance threshold) — in (0,1), typical 0.05

        - Must be in (0,1) — refuses 0, 1, negative
        - Typical 0.05 = 5% FDR
        - Smaller alpha more stringent, larger more permissive

        See JSON schema exclusiveMinimum 0 exclusiveMaximum 1.
        """,
        "advanced.dispersion_method" => """
        Dispersion estimation (NB_GLM advanced, behind Advanced Analysis)

        - parametric: fit dispersion ~ mean trend (DESeq2 default, recommended)
        - local: local regression fit (when parametric fails)
        - mean: use mean dispersion (when n small)
        - pooled: pool across genes (when n very small)
        - glmGamPoi: fast estimator from glmGamPoi package (deferred, fast)

        Heavy validation, only meaningful for NB_GLM, warning if used for other methods.

        See Love et al. 2014 and context_help('method').
        """,
        "advanced.zero_handling" => """
        Zero handling (advanced, behind Advanced Analysis, hidden unless Evidence Mode)

        - pseudocount: add pseudocount (default, safe)
        - multiplicative_replacement: multiplicative replacement (Martín-Fernández)
        - bayesian_multiplicative: Bayesian multiplicative
        - refuse: refuse to handle zeros — DANGEROUS, requires DANGER token, mathematically invalid for CLR/ILR

        See context_help('normalization.zero_policy') and Nickel ZeroHandlingContract.
        Heavy validation, refusal of meaningless, DANGER banner if refuse.
        """,
        "advanced.zero_policy" => """
        Zero policy (advanced, same as zero_handling but enum)

        See advanced.zero_handling — pseudocount, multiplicative_replacement, bayesian_multiplicative, refuse.

        Refuse is DANGEROUS and requires token, but still refused at runtime for CLR/ILR because log(0) undefined.

        See Nickel ZeroHandlingContract.
        """,
        "advanced.pseudocount" => """
        Custom pseudocount (advanced, behind Advanced Analysis, hidden unless Evidence Mode)

        - Must be >0, typical 0.5
        - <0.1 very small → extreme log-ratios for zeros — warning
        - >=1 unusual → distorts low-abundance — warning
        - Heavy validation, refusal if <=0
        - Warnings for paper writers: "pseudocount <0.1 very small will create extreme log-ratios" etc.

        See normalization.pseudocount and Nickel PseudocountContract.
        """,
        "advanced.epsilon" => """
        Epsilon for numerical stability (advanced, behind Advanced Analysis)

        - Must be in (0,1), typical 1e-6
        - >1e-3 large may affect transforms — warning
        - <1e-12 extremely small may cause underflow — warning
        - Heavy validation refusal if not in (0,1)

        See normalization.epsilon.
        """,
        "advanced.min_prevalence" => """
        Minimum prevalence filter [0,1] (advanced)

        - Feature must be present in at least this fraction of samples
        - 0.1 = present in >=10% samples — recommended to reduce multiple testing burden
        - 0 = no filter, 1 = present in 100% samples
        - Heavy validation [0,1]

        See JSON schema and Nickel PrevalenceContract.
        """,
        "advanced.min_abundance" => """
        Minimum abundance filter >=0 (advanced)

        - Feature must have at least this total abundance across samples
        - 0 = no filter
        - Heavy validation >=0
        """,
        "advanced.max_features" => """
        Maximum features to test (advanced, optional)

        - >0 or nothing
        - <=0 meaningless — refuses 0
        - >100000 excessive memory — refuses
        - <10 very small — warning will test only N features may miss biology
        - Used to limit to top N abundant/prevalent features for speed
        """,
        "advanced.min_samples_per_group" => """
        Minimum samples per group (advanced)

        - Must be >=2, recommended >=3
        - <2 refuses — need at least 2 per group for statistical test
        - <3 warning — power very low, results unreliable
        - DANGEROUS if <3 — marks config dangerous, triggers DANGER banner

        See JSON schema and Nickel.
        """,
    )

    return get(help_db, field_path, "No help available for '\$field_path'. See JSON schema and Nickel contracts. Field path examples: method, formula, metadata_columns, normalization.method, correction.method, advanced.pseudocount, advanced.epsilon, advanced.zero_policy.")
end

# --------------------------------------------------------------------------
# DANGER banner — scary for paper writers on overrides
# --------------------------------------------------------------------------

"""
    danger_banner(config)

Return a formatted warning for a configuration whose stored danger flag is set,
or `nothing` for a safe configuration. The banner lists the risky settings
recognised by this formatter and includes configuration provenance identifiers.
"""
function danger_banner(config::AnalysisConfig)
    if !is_dangerous(config)
        return nothing
    end

    reasons = String[]
    if config.correction.allow_no_correction
        push!(reasons, "BH correction disabled (method='$(config.correction.method)') — will inflate false discoveries (e.g., 1500 taxa → ~75 false positives under null at alpha=0.05)")
    end
    if config.advanced.zero_handling == "refuse" || config.advanced.zero_policy == REFUSE
        push!(reasons, "zero_handling='refuse' — will cause log(0) for CLR/ILR and biased handling for NB_GLM — mathematically invalid for CLR/ILR even with token")
    end
    if config.advanced.min_samples_per_group < 3
        push!(reasons, "min_samples_per_group=$(config.advanced.min_samples_per_group) <3 — statistical power very low, results unreliable")
    end
    if config.normalization.method == "rarefy" && config.method == NB_GLM
        push!(reasons, "rarefy + NB_GLM — rarefy discards data and NB_GLM already handles library size via size_factors — combining is questionable")
    end

    banner = """
    ╔════════════════════════════════════════════════════════════════════════════╗
    ║  ⚠️  DANGER — SCIENTIFICALLY RISKY CONFIGURATION DETECTED  ⚠️             ║
    ║  This configuration overrides safe defaults and may produce             ║
    ║  misleading or irreproducible results. Review carefully before          ║
    ║  publishing. This banner will be logged, included in provenance,       ║
    ║  and in DOI bundle.                                                    ║
    ╠════════════════════════════════════════════════════════════════════════════╣
    ║  Reasons:                                                              ║
    $(join(["║  - $r" for r in reasons], "\n"))
    ║                                                                        ║
    ║  Acknowledgment token: $(config.correction.acknowledgment_token === nothing ? config.advanced.acknowledgment_token : config.correction.acknowledgment_token) ║
    ║  Config ID: $(config.id)                                                ║
    ║  Hash: $(config.hash)                                                   ║
    ║  Method: $(METHOD_TO_STRING[config.method])                             ║
    ║  Formula: $(config.formula)                                             ║
    ║                                                                        ║
    ║  If you are writing a paper, you MUST disclose these overrides in      ║
    ║  Methods and discuss limitations. Uncorrected p-values in high-dim     ║
    ║  data are NOT publishable without strong justification.                ║
    ╚════════════════════════════════════════════════════════════════════════════╝
    """

    return banner
end

"""Return the danger flag stored on `config`."""
function is_dangerous(config::AnalysisConfig)
    return config.dangerous
end

"""
    log_danger_banner(config)

Log an informational message and return `nothing` for a safe configuration. For
a dangerous configuration, log the generated banner at error and warning levels
and return the banner text.
"""
function log_danger_banner(config::AnalysisConfig)
    banner = danger_banner(config)
    if isnothing(banner)
        @info "AnalysisConfig is safe — no DANGER banner" config_id=config.id method=METHOD_TO_STRING[config.method]
        return nothing
    else
        @error "DANGER BANNER — risky configuration" banner config_id=config.id method=METHOD_TO_STRING[config.method] hash=config.hash
        # Also log to file for paper writers
        @warn "SCARY DANGER BANNER FOR PAPER WRITERS — overrides detected, see logs and provenance" banner config_id=config.id
        return banner
    end
end

# --------------------------------------------------------------------------
# Serialization — JSON, Nickel, DEED, with standards
# --------------------------------------------------------------------------

"""Serialise `config` to the JSON representation used by the analysis API."""
function to_json(config::AnalysisConfig)
    return JSON3.write(OrderedDict(
        "schema_version" => config.schema_version,
        "id" => config.id,
        "created_at" => string(config.created_at),
        "created_by" => config.created_by,
        "method" => METHOD_TO_STRING[config.method],
        "formula" => config.formula,
        "outcome_column" => config.outcome_column,
        "metadata_columns" => config.metadata_columns,
        "normalization" => OrderedDict(
            "method" => config.normalization.method,
            "pseudocount" => config.normalization.pseudocount,
            "epsilon" => config.normalization.epsilon,
            "zero_policy" => string(config.normalization.zero_policy),
            "ilr_basis" => config.normalization.ilr_basis,
            "multiplicative_replacement_delta" => config.normalization.multiplicative_replacement_delta
        ),
        "correction" => OrderedDict(
            "method" => config.correction.method,
            "alpha" => config.correction.alpha,
            "allow_no_correction" => config.correction.allow_no_correction,
            "acknowledgment_token" => config.correction.acknowledgment_token
        ),
        "advanced" => OrderedDict(
            "dispersion_method" => config.advanced.dispersion_method,
            "zero_handling" => config.advanced.zero_handling,
            "zero_policy" => string(config.advanced.zero_policy),
            "pseudocount" => config.advanced.pseudocount,
            "epsilon" => config.advanced.epsilon,
            "min_prevalence" => config.advanced.min_prevalence,
            "min_abundance" => config.advanced.min_abundance,
            "max_features" => config.advanced.max_features,
            "min_samples_per_group" => config.advanced.min_samples_per_group,
            "robust" => config.advanced.robust,
            "acknowledgment_token" => config.advanced.acknowledgment_token
        ),
        "provenance" => config.provenance,
        "hash" => config.hash,
        "dangerous" => config.dangerous
    ))
end

"""
    from_json(json_str)

Parse a JSON analysis configuration and rebuild its nested immutable objects.
Optional fields absent from older representations receive the constructor
defaults. JSON parsing and constructor validation errors are propagated.
"""
function from_json(json_str::String)
    data = JSON3.read(json_str)

    norm_data = data.normalization
    norm = NormalizationConfig(
        method=norm_data.method,
        pseudocount=get(norm_data, :pseudocount, 0.5),
        epsilon=get(norm_data, :epsilon, 1e-6),
        zero_policy=get(norm_data, :zero_policy, "pseudocount"),
        ilr_basis=get(norm_data, :ilr_basis, nothing),
        multiplicative_replacement_delta=get(norm_data, :multiplicative_replacement_delta, nothing)
    )

    corr_data = data.correction
    corr = CorrectionConfig(
        method=corr_data.method,
        alpha=get(corr_data, :alpha, 0.05),
        allow_no_correction=get(corr_data, :allow_no_correction, false),
        acknowledgment_token=get(corr_data, :acknowledgment_token, nothing)
    )

    adv_data = data.advanced
    adv = AdvancedConfig(
        dispersion_method=get(adv_data, :dispersion_method, "parametric"),
        zero_handling=get(adv_data, :zero_handling, "pseudocount"),
        zero_policy=get(adv_data, :zero_policy, "pseudocount"),
        pseudocount=get(adv_data, :pseudocount, 0.5),
        epsilon=get(adv_data, :epsilon, 1e-6),
        min_prevalence=get(adv_data, :min_prevalence, 0.1),
        min_abundance=get(adv_data, :min_abundance, 0.0),
        max_features=get(adv_data, :max_features, nothing),
        min_samples_per_group=get(adv_data, :min_samples_per_group, 3),
        robust=get(adv_data, :robust, false),
        acknowledgment_token=get(adv_data, :acknowledgment_token, nothing)
    )

    prov = OrderedDict{String,Any}()
    if haskey(data, :provenance)
        for (k,v) in data.provenance
            prov[string(k)] = v
        end
    end

    cfg = AnalysisConfig(
        schema_version=get(data, :schema_version, SCHEMA_VERSION),
        id=data.id,
        created_at=DateTime(data.created_at),
        created_by=get(data, :created_by, "anonymous"),
        method=data.method,
        formula=data.formula,
        outcome_column=get(data, :outcome_column, nothing),
        metadata_columns=Vector{String}(data.metadata_columns),
        normalization=norm,
        correction=corr,
        advanced=adv,
        provenance=prov,
        hash=get(data, :hash, nothing),
        dangerous=get(data, :dangerous, nothing)
    )

    return cfg
end

"""Render `config` as the module's generated Nickel source text."""
function to_nickel(config::AnalysisConfig)
    # Nickel contract from hyperpolymath/standards 1-formats/k9/*.ncl style
    # Uses TagOrString for enums, contracts for validation
    return """
    # SPDX-License-Identifier: AGPL-3.0-only
    # AnalysisConfig Nickel — $(config.id) — generated via AnalysisConfig.to_nickel()
    # Schema version $(config.schema_version), hash $(config.hash)
    # From hyperpolymath/standards: 1-formats/k9/*.ncl and .machine_readable/contractiles/_base.ncl

    let DANGER_TOKEN = "$(DANGER_ACK_TOKEN)" in

    {
      schema_version = "$(config.schema_version)",
      id = "$(config.id)",
      created_at = "$(config.created_at)",
      created_by = "$(config.created_by)",
      method = '$(METHOD_TO_STRING[config.method])',
      formula = "$(config.formula)" | ValidFormula,
      outcome_column = $(isnothing(config.outcome_column) ? "null" : "\"$(config.outcome_column)\""),
      metadata_columns = [$(join(["\"\$c\"" for c in config.metadata_columns], ", "))],

      normalization = {
        method = '$(config.normalization.method)',
        pseudocount = $(config.normalization.pseudocount) | PseudocountContract,
        epsilon = $(config.normalization.epsilon),
        zero_policy = '$(string(config.normalization.zero_policy))',
        ilr_basis = $(isnothing(config.normalization.ilr_basis) ? "null" : "'$(config.normalization.ilr_basis)'"),
      } | MethodNormalizationCompatibility,

      correction = {
        method = '$(config.correction.method)',
        alpha = $(config.correction.alpha),
        allow_no_correction = $(config.correction.allow_no_correction ? "true" : "false"),
        acknowledgment_token = $(isnothing(config.correction.acknowledgment_token) ? "null" : "\"$(config.correction.acknowledgment_token)\""),
      } | CorrectionContract,

      advanced = {
        dispersion_method = '$(config.advanced.dispersion_method)',
        zero_handling = '$(config.advanced.zero_handling)',
        pseudocount = $(config.advanced.pseudocount) | PseudocountContract,
        epsilon = $(config.advanced.epsilon),
        min_prevalence = $(config.advanced.min_prevalence) | PrevalenceContract,
        min_abundance = $(config.advanced.min_abundance),
        max_features = $(isnothing(config.advanced.max_features) ? "null" : string(config.advanced.max_features)),
        min_samples_per_group = $(config.advanced.min_samples_per_group),
        robust = $(config.advanced.robust ? "true" : "false"),
      } | ZeroHandlingContract,

      provenance = {
        hash = "$(config.hash)",
        dangerous = $(config.dangerous ? "true" : "false"),
      },

      hash = "$(config.hash)",
      dangerous = $(config.dangerous ? "true" : "false"),
    }
    """
end

"""
    from_nickel(nickel_str)

Build a minimal configuration from the first method and formula assignments in
Nickel-like source. This is a partial parser: metadata and most settings use
defaults. Raise `ArgumentError` when either assignment cannot be found or the
derived configuration is invalid.
"""
function from_nickel(nickel_str::String)
    # For now, parse via simple regex — full Nickel evaluation would require nickel binary
    # This is a placeholder that extracts method and formula and validates via our validators
    # Real implementation would call `nickel eval` via pipeline tools
    m = match(r"method\s*=\s*'(\w+)'", nickel_str)
    f = match(r"formula\s*=\s*\"([^\"]+)\"", nickel_str)
    if isnothing(m) || isnothing(f)
        throw(ArgumentError("Failed to parse Nickel: could not find method and formula — see to_nickel() output"))
    end
    method_str = m.captures[1]
    formula_str = f.captures[1]
    # For demo, create minimal config — real would parse full Nickel
    return AnalysisConfig(
        method=method_str,
        formula=formula_str,
        metadata_columns=["group"], # minimal, would be parsed from Nickel in real
        normalization=NormalizationConfig(method=method_str == "clr_lm" ? "clr" : method_str == "ilr_lm" ? "ilr" : "size_factors"),
    )
end

"""
    validate_nickel(nickel_str)

Return messages for required field or contract markers missing from generated
Nickel text. This performs textual checks only; it does not parse or evaluate
Nickel.
"""
function validate_nickel(nickel_str::String)
    # Validate via contracts — for now check that it contains required fields and contracts pass
    errors = String[]
    if !occursin("schema_version", nickel_str)
        push!(errors, "Nickel missing schema_version")
    end
    if !occursin("method", nickel_str)
        push!(errors, "Nickel missing method")
    end
    if !occursin("formula", nickel_str)
        push!(errors, "Nickel missing formula")
    end
    if !occursin("ValidFormula", nickel_str)
        push!(errors, "Nickel missing ValidFormula contract — see hyperpolymath/standards 1-formats/k9/")
    end
    if !occursin("CorrectionContract", nickel_str)
        push!(errors, "Nickel missing CorrectionContract (BH mandatory)")
    end
    return errors
end

"""Render `config` as a `repo-deed` document using DEED syntax."""
function to_deed(config::AnalysisConfig)
    # DEED from hyperpolymath/standards 1-formats/deed/spec/DEED-GRAMMAR-SPEC.adoc v0.2.0
    # :schema-version first, only () brackets, #t/#f booleans, :kebab-case keywords, #u5 UUID5, SPDX header mandatory
    dangerous_bool = config.dangerous ? "#t" : "#f"
    allow_no_corr_bool = config.correction.allow_no_correction ? "#t" : "#f"
    robust_bool = config.advanced.robust ? "#t" : "#f"

    return """
    ;; SPDX-License-Identifier: AGPL-3.0-only
    ;; SPDX-FileCopyrightText: 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
    ;; AnalysisConfig DEED — $(config.id) — generated via AnalysisConfig.to_deed()
    ;; Schema version $(config.schema_version), hash $(config.hash)
    ;; From hyperpolymath/standards: 1-formats/deed/spec/DEED-GRAMMAR-SPEC.adoc v0.2.0
    ;; :schema-version first, only () brackets, #t/#f booleans, :kebab-case keywords, #u5 UUID5

    (repo-deed
      :schema-version "$(config.schema_version)"
      :canonical-name "analysis-config-$(config.id)"
      :beholding-chora #u5"estate/chora"

      (method
        :name "$(METHOD_TO_STRING[config.method])"
        :formula "$(config.formula)"
        :outcome-column "$(isnothing(config.outcome_column) ? "" : config.outcome_column)"
        :metadata-columns ($(join(["\"\$c\"" for c in config.metadata_columns], " "))))

      (normalization
        :method "$(config.normalization.method)"
        :pseudocount $(config.normalization.pseudocount)
        :epsilon $(config.normalization.epsilon)
        :zero-policy "$(string(config.normalization.zero_policy))"
        :ilr-basis "$(isnothing(config.normalization.ilr_basis) ? "" : config.normalization.ilr_basis)"
        :multiplicative-replacement-delta $(isnothing(config.normalization.multiplicative_replacement_delta) ? "0" : string(config.normalization.multiplicative_replacement_delta)))

      (correction
        :method "$(config.correction.method)"
        :alpha $(config.correction.alpha)
        :allow-no-correction $(allow_no_corr_bool)
        :acknowledgment-token "$(isnothing(config.correction.acknowledgment_token) ? "" : config.correction.acknowledgment_token)")

      (advanced
        :dispersion-method "$(config.advanced.dispersion_method)"
        :zero-handling "$(config.advanced.zero_handling)"
        :zero-policy "$(string(config.advanced.zero_policy))"
        :pseudocount $(config.advanced.pseudocount)
        :epsilon $(config.advanced.epsilon)
        :min-prevalence $(config.advanced.min_prevalence)
        :min-abundance $(config.advanced.min_abundance)
        :max-features $(isnothing(config.advanced.max_features) ? "0" : string(config.advanced.max_features))
        :min-samples-per-group $(config.advanced.min_samples_per_group)
        :robust $(robust_bool)
        :acknowledgment-token "$(isnothing(config.advanced.acknowledgment_token) ? "" : config.advanced.acknowledgment_token)")

      (provenance
        :id "$(config.id)"
        :hash "$(config.hash)"
        :created-at "$(config.created_at)"
        :created-by "$(config.created_by)"
        :dangerous $(dangerous_bool)
        :schema-version "$(config.schema_version)")

      (warrant
        :evidence-type "AnalysisConfig"
        :soundness "$(config.correction.allow_no_correction ? "UNSOUND — BH disabled" : "BH mandatory, requires acknowledgment token for override")"
        :fiber "Echo of raw counts through $(config.normalization.method) transform"
        :epistemic-status "$(config.dangerous ? "present_in_some_admissible_world" : "present_in_every_admissible_world")")

      (doi-bundle
        :title "MetaManifold Analysis Bundle $(config.id)"
        :license "CC-BY-4.0"
        :authors ("$(config.created_by)")
        :description "Differential abundance $(METHOD_TO_STRING[config.method]) formula $(config.formula) BH $(config.correction.method) DOI-ready")

      (context-help
        :method "Analysis method explicit no auto-selection v1 nb_glm clr_lm ilr_lm logistic"
        :formula "R-style formula e.g. ~ group must reference only metadata_columns forbids ; backtick dollar"
        :correction "BH mandatory in v1 any override triggers DANGER banner requires acknowledgment token $(DANGER_ACK_TOKEN)"
        :normalization "Normalization must be compatible with method nb_glm allows none/rarefy/size_factors/relative/TSS/CLR/ILR"
        :advanced "All advanced options behind Advanced Analysis expander hidden unless Evidence Mode heavy validation refusal meaningless"))
    """
end

"""
    validate_deed(deed_str)

Return messages for required DEED markers, forbidden bracket characters, and a
missing DEED boolean. These are textual checks rather than full DEED parsing.
"""
function validate_deed(deed_str::String)
    errors = String[]
    if !occursin(":schema-version", deed_str)
        push!(errors, "DEED missing :schema-version (must be first per DEED-GRAMMAR-SPEC v0.2.0)")
    end
    if !occursin("repo-deed", deed_str)
        push!(errors, "DEED missing repo-deed head (filename dispatch *_chora.deed → repo-deed)")
    end
    if !occursin("SPDX-License-Identifier:", deed_str)
        push!(errors, "DEED missing SPDX-License-Identifier: header (mandatory per spec)")
    end
    if occursin(r"[\[\]\{\}]", deed_str)
        # DEED allows only () brackets per spec
        # But our template uses only (), so check for [] {} which are forbidden
        # However JSON arrays inside strings are okay, but brackets outside parens are forbidden
        # Simple check: if [] or {} appear outside of quoted strings, it's invalid
        # For now, just warn if [] {} appear at all (since our template uses () only)
        # Actually our template uses []? No, uses () only, so [] {} would be invalid
        # We'll check for [ or ] or { or } not inside quotes — simplified
        push!(errors, "DEED contains forbidden brackets [] or {} — only () allowed per DEED-GRAMMAR-SPEC v0.2.0")
    end
    if !occursin("#t", deed_str) && !occursin("#f", deed_str)
        push!(errors, "DEED missing #t/#f booleans (should use #t/#f per spec, not true/false)")
    end
    return errors
end

# --------------------------------------------------------------------------
# DOI-ready JSON manifest bundles — DataCite
# --------------------------------------------------------------------------

"""
    create_doi_bundle(config, [result]; output_dir, authors, title, license, description)

Create or update `output_dir` with DataCite metadata, JSON, Nickel, DEED,
provenance, and content-hash files, then return the directory path. An optional
result is embedded in the DataCite document. Dangerous configurations also write
`DANGER_BANNER.txt` and emit a warning; existing files with the same names are
overwritten.
"""
function create_doi_bundle(config::AnalysisConfig, result::Union{AnalysisResult,Nothing}=nothing; output_dir::String="doi_bundle_$(config.id)", authors::Vector{String}=String[], title::String="MetaManifold Analysis Bundle", license::String="CC-BY-4.0", description::String="Differential abundance analysis")
    mkpath(output_dir)

    # DataCite JSON
    datacite = OrderedDict{String,Any}(
        "id" => config.id,
        "type" => "Dataset",
        "titles" => [OrderedDict("title" => title)],
        "creators" => [OrderedDict("name" => a) for a in (isempty(authors) ? [config.created_by] : authors)],
        "descriptions" => [OrderedDict("description" => description, "descriptionType" => "Abstract")],
        "publicationYear" => year(config.created_at),
        "publisher" => "MetaManifold-WebUI",
        "resourceType" => OrderedDict("resourceTypeGeneral" => "Dataset", "resourceType" => "AnalysisConfig"),
        "subjects" => [
            OrderedDict("subject" => "microbiome"),
            OrderedDict("subject" => "differential abundance"),
            OrderedDict("subject" => METHOD_TO_STRING[config.method]),
            OrderedDict("subject" => "BH correction"),
        ],
        "formats" => ["application/json", "text/nickel", "text/deed"],
        "version" => config.schema_version,
        "rightsList" => [OrderedDict("rights" => license)],
        "dates" => [OrderedDict("date" => string(config.created_at), "dateType" => "Created")],
        "relatedIdentifiers" => [
            OrderedDict("relatedIdentifier" => config.hash, "relatedIdentifierType" => "SHA256", "relationType" => "IsIdenticalTo"),
        ],
        "schemaVersion" => "http://datacite.org/schema/kernel-4",
        "config" => JSON3.read(to_json(config)),
        "provenance" => config.provenance,
        "dangerous" => config.dangerous,
        "warrant" => OrderedDict(
            "evidence_type" => "AnalysisConfig",
            "soundness" => config.correction.allow_no_correction ? "UNSOUND — BH disabled" : "BH mandatory",
            "fiber" => "Echo of raw counts through $(config.normalization.method) transform",
            "epistemic_status" => config.dangerous ? "present_in_some_admissible_world" : "present_in_every_admissible_world"
        )
    )

    if !isnothing(result)
        datacite["result"] = OrderedDict(
            "id" => result.id,
            "config_id" => result.config_id,
            "config_hash" => result.config_hash,
            "hash" => result.hash,
            "method" => METHOD_TO_STRING[result.method],
            "results" => result.results
        )
        datacite["relatedIdentifiers"] = vcat(datacite["relatedIdentifiers"], [
            OrderedDict("relatedIdentifier" => result.hash, "relatedIdentifierType" => "SHA256", "relationType" => "HasPart")
        ])
    end

    # Write files
    open(joinpath(output_dir, "datacite.json"), "w") do io
        JSON3.write(io, datacite)
    end

    open(joinpath(output_dir, "analysis_config.json"), "w") do io
        write(io, to_json(config))
    end

    open(joinpath(output_dir, "analysis_config.ncl"), "w") do io
        write(io, to_nickel(config))
    end

    open(joinpath(output_dir, "analysis_config_chora.deed"), "w") do io
        write(io, to_deed(config))
    end

    # Provenance file
    open(joinpath(output_dir, "provenance.json"), "w") do io
        JSON3.write(io, config.provenance)
    end

    # DANGER banner if dangerous
    if config.dangerous
        banner = danger_banner(config)
        open(joinpath(output_dir, "DANGER_BANNER.txt"), "w") do io
            write(io, banner)
        end
        @warn "DOI bundle contains DANGER banner — config is dangerous" output_dir config_id=config.id
    end

    # Content-addressed hash file
    open(joinpath(output_dir, "content_hash.txt"), "w") do io
        write(io, config.hash)
    end

    @info "Created DOI-ready bundle" output_dir config_id=config.id hash=config.hash dangerous=config.dangerous

    return output_dir
end

# --------------------------------------------------------------------------
# Epistemic bridge — present_in_every_admissible_world (finite model)
# --------------------------------------------------------------------------

"""
    present_in_every_admissible_world(candidates, query)

Return whether `query` is true for every candidate world. An empty candidate
collection therefore returns `true`, following `all` semantics.
"""
function present_in_every_admissible_world(candidates::Vector, query::Function)
    # Finite model from residual-evidence-types: Case inhabited, Holds c P = (x:Candidate) -> P (fst x)
    # Returns true if query holds for every candidate world consistent with observation
    # Example: r = u + n = 2 with NoiseBound n ≤ 1 gives presence without identification (u ≠ 0 but u=1 or 2)
    return all(c -> query(c), candidates)
end

end # module AnalysisConfigImpl
