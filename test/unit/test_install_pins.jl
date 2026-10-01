using Test
using YAML
using JSON3
using SHA

## Pinning of external tools and runtimes
#
# These tests police the pins themselves rather than any behaviour of the pipeline.
# They exist because the failure they guard against is silent: an installer that
# takes whatever upstream published most recently still installs, still passes, and
# still produces results, and nothing anywhere records that two runs were made by
# different software. The only defence is to assert that the pins are declared, that
# they are declared in one place, and that every copy of them agrees.

const REPO_ROOT   = normpath(joinpath(@__DIR__, "..", ".."))
const PINS_PATH   = joinpath(REPO_ROOT, "config", "defaults", "tool_versions.yml")
const INSTALL_JL  = joinpath(REPO_ROOT, "install.jl")
const CI_PATH     = joinpath(REPO_ROOT, ".github", "workflows", "ci.yml")
const MANIFEST    = joinpath(REPO_ROOT, "Manifest.toml")
const RENV_LOCK   = joinpath(REPO_ROOT, "renv.lock")

# The tools the pipeline shells out to. A tool absent from the pin file is a tool
# nobody has decided the version of.
const EXPECTED_TOOLS = ["cutadapt", "fastqc", "multiqc", "vsearch", "cd_hit_est", "swarm",
                        "mafft", "trimal", "iqtree", "raxml", "gappa"]

# The platforms install.jl can resolve an archive for, as "<OS_TYPE>-<ARCH_STR>".
const EXPECTED_PLATFORMS = ["linux-x86_64", "linux-aarch64", "macos-x86_64", "macos-aarch64"]

is_sha256(s) = s isa AbstractString && occursin(r"^[0-9a-f]{64}$", s)

"""Return the first step of `job` whose `uses:` names `action`, or nothing."""
function ci_step_using(job, action)
    for step in get(job, "steps", [])
        startswith(get(step, "uses", ""), action) && return step
    end
    return nothing
end


@testset "install pins" begin
    @test isfile(PINS_PATH)
    pins = YAML.load_file(PINS_PATH)

    @testset "every tool is pinned" begin
        @test pins["schema_version"] == 1
        @test sort(collect(keys(pins["tools"]))) == sort(EXPECTED_TOOLS)

        for tool in EXPECTED_TOOLS
            rec = pins["tools"][tool]
            version = get(rec, "version", nothing)
            @test version isa AbstractString
            @test !isempty(version)
            # "latest" is the absence of a pin wearing the costume of one.
            @test !occursin("latest", lowercase(version))
        end
    end

    @testset "every archive is checksummed" begin
        for tool in EXPECTED_TOOLS
            rec = pins["tools"][tool]
            rec["source"] == "pypi" && continue

            archives = rec["archives"]
            @test !isempty(archives)
            for (platform, archive) in archives
                url = get(archive, "url", nothing)
                @test url isa AbstractString
                @test startswith(url, "https://")
                # The version must be legible in the URL, or the pin and the artefact
                # it names can part company without anyone noticing.
                @test occursin(rec["version"], url)

                # A null checksum is permitted by the file's own documentation, but only
                # as an explicit, visible admission. It must never pass unremarked.
                sha = get(archive, "sha256", nothing)
                @test is_sha256(sha)
                sha === nothing && @warn "No SHA256 pinned for $tool on $platform; " *
                                         "that download cannot be verified."
            end
        end
    end

    @testset "binary tools are pinned for every supported platform" begin
        # FastQC is Java and cd-hit is built from source, so both are platform-neutral;
        # vsearch and swarm ship per-platform binaries and must cover what install.jl
        # will ask for.
        for tool in ("vsearch", "swarm", "trimal", "iqtree")
            @test sort(collect(keys(pins["tools"][tool]["archives"]))) == sort(EXPECTED_PLATFORMS)
        end
        for tool in ("fastqc", "cd_hit_est", "mafft", "raxml", "gappa")
            @test collect(keys(pins["tools"][tool]["archives"])) == ["any"]
        end
    end

    @testset "install.jl holds no versions of its own" begin
        src = read(INSTALL_JL, String)

        # The regression that prompted all of this: vsearch and swarm were installed
        # from whatever GitHub called newest that day.
        @test !occursin("releases/latest", src)
        @test !occursin("api.github.com", src)

        # Download URLs live in the pin file, so install.jl must carry none. Were one to
        # creep back, the pin file would cease to be the single source of truth and the
        # two would drift apart in silence.
        @test !occursin("https://github.com/torognes", src)
        @test !occursin("bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v", src)
        @test !occursin("github.com/weizhongli", src)
        @test !occursin("mafft.cbrc.jp", src)
        @test !occursin("github.com/iqtree", src)

        @test occursin("tool_versions.yml", src)
    end

    @testset "runtime pins agree with the lockfiles" begin
        # The pin file does not lead the lockfiles; it restates them so that CI and the
        # installer can read a version the lockfiles hold in formats they cannot parse.
        julia_pin = pins["runtimes"]["julia"]["version"]
        manifest  = read(MANIFEST, String)
        m = match(r"julia_version\s*=\s*\"([^\"]+)\"", manifest)
        @test m !== nothing
        @test julia_pin == m[1]

        renv = JSON3.read(read(RENV_LOCK, String))
        @test pins["runtimes"]["r"]["version"] == renv["R"]["Version"]
        @test pins["runtimes"]["r"]["bioconductor"] == renv["Bioconductor"]["Version"]

        # The apt revision is how R is actually pinned on the CI runner, so it must name
        # the same R the lockfile does.
        @test startswith(pins["runtimes"]["r"]["apt_version"], pins["runtimes"]["r"]["version"] * "-")
    end

    @testset "Julia is bootstrapped from a pinned, checksummed archive" begin
        julia = pins["runtimes"]["julia"]
        @test sort(collect(keys(julia["archives"]))) == sort(EXPECTED_PLATFORMS)
        for (platform, archive) in julia["archives"]
            @test startswith(archive["url"], "https://")
            @test occursin(julia["version"], archive["url"])
            @test is_sha256(archive["sha256"])
        end

        # install.sh used to pipe a live installer into sh: no version, no checksum.
        src = read(joinpath(REPO_ROOT, "install.sh"), String)
        @test !occursin("install.julialang.org", src)
        @test !occursin(r"\|\s*sh\b", src)
        @test occursin("sha256", src)
    end

    @testset "the installer refuses what it cannot verify" begin
        # install.jl guards its own main(), so loading it here resolves the pins and
        # defines the helpers without installing anything.
        Installer = Module(:Installer)
        Base.include(Installer, INSTALL_JL)

        @test Installer.pinned_version("vsearch") == pins["tools"]["vsearch"]["version"]

        # Platform-specific archives are selected by platform; platform-neutral ones fall
        # back to the "any" entry rather than failing.
        vsearch = Installer.pinned_archive("vsearch")
        @test vsearch["url"] == pins["tools"]["vsearch"]["archives"][Installer.PLATFORM]["url"]
        @test Installer.pinned_archive("fastqc")["url"] ==
              pins["tools"]["fastqc"]["archives"]["any"]["url"]

        mktempdir() do dir
            payload = joinpath(dir, "payload.tar.gz")
            write(payload, "the bytes an installer was promised")
            actual = bytes2hex(open(SHA.sha256, payload))

            # The honest case: the hash agrees, and the file survives.
            @test Installer.verify_sha256(payload, actual, "payload.tar.gz") === nothing
            @test isfile(payload)

            # The dishonest case. Refusal must be loud, and the suspect file must not be
            # left lying about for the next run to pick up and trust.
            wrong = "0" ^ 64
            @test_throws ErrorException Installer.verify_sha256(payload, wrong, "payload.tar.gz")
            @test !isfile(payload)

            # An absent checksum is not permission to install unverified bytes.
            write(payload, "anything at all")
            @test_throws ErrorException Installer.verify_sha256(payload, nothing, "payload.tar.gz")
        end
    end

    @testset "CI installs what is pinned" begin
        ci = YAML.load_file(CI_PATH)
        job = ci["jobs"]["test"]

        # Julia is the one version CI cannot read from the pin file, because nothing can
        # be read before Julia exists. It is therefore duplicated, and this is the test
        # that makes the duplication safe.
        #
        # Read from the `Set up Julia` step rather than from a matrix: the matrix was
        # removed because GitHub appends a matrix combination to the posted check name
        # (see "a required check name is a stable identifier" below). The step is located
        # by its `uses:` rather than by index, so reordering the steps cannot make this
        # assertion quietly vanish.
        setup = ci_step_using(job, "julia-actions/setup-julia")
        @test setup !== nothing
        @test setup["with"]["version"] == pins["runtimes"]["julia"]["version"]

        # A floating runner would carry the R apt pin, which names a 24.04 build, off to
        # whatever the next LTS ships.
        @test job["runs-on"] == "ubuntu-24.04"
        @test occursin("2404", pins["runtimes"]["r"]["apt_version"])

        # Everything else CI installs must be read from the pin file at run time rather
        # than written out beside it.
        steps = ci["jobs"]["test"]["steps"]
        runs  = join([get(step, "run", "") for step in steps], "\n")
        @test occursin("tool_versions.yml", runs)
        @test occursin("sha256sum -c", runs)

        for tool in ("vsearch", "swarm", "fastqc", "multiqc")
            @test !occursin(pins["tools"][tool]["version"], runs)
        end

        # Every tool the pipeline calls is installed. The check looks for each pin's
        # variable where a step uses it, because the pin-export step names every tool.
        for (tool, ref) in (("cutadapt", raw"$CUTADAPT_SPEC"),
                            ("multiqc",  raw"$MULTIQC_SPEC"),
                            ("fastqc",   raw"$FASTQC_URL"),
                            ("vsearch",  raw"$VSEARCH_URL"),
                            ("swarm",    raw"$SWARM_URL"))
            @test occursin(ref, runs)
        end

        # config/ci/tools.yml resolves every tool by bare name.
        for tool in ("vsearch", "swarm", "fastqc")
            @test occursin("/usr/local/bin/$tool", runs)
        end

        # cd-hit is the one tool CI installs from apt without a pin.
        @test occursin("apt-get install -y cd-hit", runs)

        # Archive downloads go through fetch_pinned, which retries TLS failures and
        # names the cause when a download still fails. Comments are stripped so the
        # checks read the code.
        helper = joinpath(REPO_ROOT, "scripts", "ci", "fetch_pinned.sh")
        @test isfile(helper)
        fetch_code = join([line for line in eachline(helper)
                           if !startswith(strip(line), "#")], "\n")
        @test occursin("--retry-all-errors", fetch_code)
        @test occursin("::error::TLS verification failed", fetch_code)
        @test occursin("tool_versions.yml", fetch_code)
        @test !occursin("|| true", fetch_code)
        @test !occursin("set +e", fetch_code)

        @test occursin("source scripts/ci/fetch_pinned.sh", runs)
        for ref in (raw"$VSEARCH_URL", raw"$SWARM_URL", raw"$FASTQC_URL")
            @test occursin("fetch_pinned \"$ref\"", runs)
        end
    end

    @testset "a required check name is a stable identifier" begin
        ci = YAML.load_file(CI_PATH)

        # A branch ruleset matches a required status check by the DISPLAY NAME GitHub
        # posts for the job. A renamed check does not fail -- it is simply absent, and a
        # required check that never reports can never be satisfied, so every pull request
        # deadlocks until an admin bypasses the rule. The bump that triggers it is a
        # one-line edit to this very file's subject matter.
        #
        # TWO conditions are needed for GitHub to post the `name:` field verbatim, and
        # the first without the second is the trap this file exists to close:
        #
        #   1. the name must not interpolate a pin; and
        #   2. the job must not have a `strategy.matrix`.
        #
        # MEASURED 2026-09-21 on hyperpolymath/MetaManifold-WebUI#39: a job named exactly `Julia tests` with a 1x1
        # matrix posted `Julia tests (1.12.5, ubuntu-24.04)`. GitHub appends the matrix
        # combination whenever the name does not already reference the matrix, so a 1x1
        # matrix is still a matrix and the version was still embedded. A guard asserting
        # only (1) PASSED on that broken fix, because it asked about the YAML field while
        # the ruleset reads the rendered check name.
        #
        # Asserted over EVERY job rather than the one job we happen to require today:
        # the rule is "a job name is a stable identifier", and a rule enforced at each
        # door in turn is a rule that a new door escapes. If a matrix is ever genuinely
        # needed, the required context must move to a matrix-free aggregator job and this
        # assertion be re-scoped to that job -- not deleted.
        for (id, job) in ci["jobs"]
            name = get(job, "name", id)
            @test !occursin("\${{", name)
            @test !haskey(get(job, "strategy", Dict{String,Any}()), "matrix")
        end
    end
end
