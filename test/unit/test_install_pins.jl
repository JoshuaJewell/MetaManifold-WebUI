# SPDX-License-Identifier: AGPL-3.0-only
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
const EXPECTED_TOOLS = ["cutadapt", "fastqc", "multiqc", "vsearch", "cd_hit_est", "swarm"]

# The platforms install.jl can resolve an archive for, as "<OS_TYPE>-<ARCH_STR>".
const EXPECTED_PLATFORMS = ["linux-x86_64", "linux-aarch64", "macos-x86_64", "macos-aarch64"]

is_sha256(s) = s isa AbstractString && occursin(r"^[0-9a-f]{64}$", s)

"""Return the first step of `job` whose `name:` is `name`, or nothing."""
function ci_step_named(job, name)
    for step in get(job, "steps", [])
        get(step, "name", nothing) == name && return step
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
        for tool in ("vsearch", "swarm")
            @test sort(collect(keys(pins["tools"][tool]["archives"]))) == sort(EXPECTED_PLATFORMS)
        end
        for tool in ("fastqc", "cd_hit_est")
            @test collect(keys(pins["tools"][tool]["archives"])) == ["any"]
        end
    end

    @testset "the FastQC archive is vendored, and stays vendored" begin
        # Guarded as a decision rather than described as a fact, because the failure
        # mode is silent: a URL that quietly points back at a third party still
        # downloads, still verifies against the checksum, and gives no sign that CI is
        # once again dependent on someone else's certificate renewal.
        #
        # MEASURED 2026-09-22: CI went red on a commit that touched nothing near FastQC
        # because www.bioinformatics.babraham.ac.uk served an expired TLS certificate.
        # FastQC publishes no release assets upstream -- all six of its GitHub releases
        # carry none -- so the pin had nowhere else to point. The archive is now served
        # from this repository's own releases (see docs/compliance/vendored-archives.md):
        # the same bytes, kept against the checksum that was already pinned, so only the
        # serving host changed.
        #
        # Repointing this at a third-party host again is a decision somebody should make
        # deliberately, so it fails here by name instead of drifting back in silence.
        fastqc_url = pins["tools"]["fastqc"]["archives"]["any"]["url"]
        @test startswith(fastqc_url,
                         "https://github.com/hyperpolymath/MetaManifold-WebUI/releases/download/")
        @test !occursin("babraham.ac.uk", fastqc_url)

        # The vendored copy is only defensible while it is provably the upstream
        # artefact. This is that checksum, unchanged from the original upstream pin.
        @test pins["tools"]["fastqc"]["archives"]["any"]["sha256"] ==
              "5f4dba8780231a25a6b8e11ab2c238601920c9704caa5458d9de559575d58aa7"
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
        # by its name rather than by index, so reordering the steps cannot make this
        # assertion quietly vanish.
        #
        # It used to be located by its `uses: julia-actions/setup-julia`, but that
        # action is refused by the repository's Actions allow-list (see "every GitHub
        # Action satisfies the repository Actions policy" below), so Julia is now
        # installed by a pinned `run:` step whose JULIA_VERSION carries the copy this
        # testset exists to keep honest.
        setup = ci_step_named(job, "Set up Julia")
        @test setup !== nothing
        setup_run = setup === nothing ? "" : get(setup, "run", "")
        julia_pin = match(r"JULIA_VERSION=\"([^\"]+)\"", setup_run)
        @test julia_pin !== nothing
        if julia_pin !== nothing
            @test julia_pin[1] == pins["runtimes"]["julia"]["version"]
        end

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

        # fastqc and multiqc joined this list with issue #30. Leaving them out would
        # have made the guard vacuous for exactly the two tools it had just been
        # extended to cover -- the shape that let the gap open in the first place.
        for tool in ("vsearch", "swarm", "fastqc", "multiqc")
            @test !occursin(pins["tools"][tool]["version"], runs)
        end

        # And the gap itself: every tool the pipeline shells out to must actually be
        # installed by the workflow, which is what issue #30 found missing for fastqc
        # and multiqc.
        #
        # The assertion is on the DEREFERENCE of each exported pin, not on the tool's
        # name. Written the obvious way -- does "fastqc" appear anywhere in the run
        # blocks -- this guard passes with the install step deleted, because the
        # pin-export block above mentions every tool by name whether or not anything
        # installs it. That version was written here, and deleting the fastqc step did
        # not redden it. `$FASTQC_URL` appears only where a step consumes the pin.
        # raw"..." rather than "\$...": the string is a shell variable reference to be
        # found verbatim in the YAML, never a Julia interpolation, and config/ci/lint_source.jl
        # rejects an escaped `\$` in an interpolating string precisely to force that
        # distinction to be stated. The two forms carry the same bytes; only one says why.
        for (tool, ref) in (("cutadapt", raw"$CUTADAPT_SPEC"),
                            ("multiqc",  raw"$MULTIQC_SPEC"),
                            ("fastqc",   raw"$FASTQC_URL"),
                            ("vsearch",  raw"$VSEARCH_URL"),
                            ("swarm",    raw"$SWARM_URL"))
            @test occursin(ref, runs)
        end

        # The two archive tools that are not a single binary still have to end up on
        # PATH, because config/ci/tools.yml resolves every tool by bare name.
        for tool in ("vsearch", "swarm", "fastqc")
            @test occursin("/usr/local/bin/$tool", runs)
        end

        # cd-hit is the one tool CI does not pin, so it is matched on the apt line.
        @test occursin("apt-get install -y cd-hit", runs)

        # The archive downloads retry, and say what failed when they still do not.
        #
        # MEASURED 2026-09-22T07:09Z: main went red on commit 4df7881 -- a commit that
        # touched nothing in this area -- because the FastQC host's TLS certificate had
        # expired. The log read "curl: (60) SSL certificate problem: certificate has
        # expired" and then "##[error]Process completed with exit code 60", which reads
        # like a code failure until someone opens the log and knows what 60 means.
        #
        # curl does not retry a certificate failure unless --retry-all-errors is given,
        # so `--retry` alone would not have covered the class that bit. And skipping the
        # tool is not available either: issue #30 exists precisely because fastqc was
        # silently absent. So the installs share one helper that retries the transient
        # class and annotates the cause when the failure is real.
        helper = joinpath(REPO_ROOT, "scripts", "ci", "fetch_pinned.sh")
        @test isfile(helper)

        # Comments stripped before asserting, and that is not incidental. The obvious
        # version of this guard -- does "--retry-all-errors" appear in the file --
        # PASSES with the flag deleted from the command, because the paragraph
        # explaining the flag names it. MUTATION-TESTED 2026-09-22: removing the flag
        # from the curl invocation left this testset green until the assertion read the
        # code alone. The same trap, in the same file, for the same reason as the
        # `$FASTQC_URL` note above: a guard that can be satisfied by prose guards prose.
        fetch_code = join([line for line in eachline(helper)
                           if !startswith(strip(line), "#")], "\n")

        # The load-bearing flag: this is the fix, not a detail.
        @test occursin("--retry-all-errors", fetch_code)
        # A diagnosis rather than a bare exit code, and a lever to pull: the next
        # reader is told it is the certificate and where the pin lives.
        @test occursin("::error::TLS verification failed", fetch_code)
        @test occursin("tool_versions.yml", fetch_code)
        # It reports and then fails; it never converts a failure into a pass.
        @test !occursin("|| true", fetch_code)
        @test !occursin("set +e", fetch_code)

        # Every pinned archive goes through it, so a tool added later cannot arrive
        # with a bare curl that dies on the first TLS hiccup.
        @test occursin("source scripts/ci/fetch_pinned.sh", runs)
        for ref in (raw"$VSEARCH_URL", raw"$SWARM_URL", raw"$FASTQC_URL")
            @test occursin("fetch_pinned \"$ref\"", runs)
        end
    end

    @testset "every GitHub Action satisfies the repository Actions policy" begin
        # MEASURED 2026-09-25 22:19 UTC through 2026-09-27: every workflow run
        # failed with `startup_failure` and zero jobs because the repository
        # began enforcing an Actions allow-list. The run annotation reads:
        # actions must be "from a repository owned by hyperpolymath, created by
        # GitHub, verified in the GitHub Marketplace, or match the pattern"
        # configured for the repository, and "all actions must also be pinned
        # to a full-length commit SHA". A violation is not a red job -- it is a
        # job that never starts, so nothing INSIDE the workflow can report it.
        # This testset is where that policy becomes visible to `Pkg.test`.
        #
        # Local actions (`./...`) are part of the checked-out tree: always
        # allowed, and there is no remote ref to pin.
        github_owned_owners = ("actions", "github", "hyperpolymath")
        # Marketplace-verified third parties used by this repository. Extend
        # deliberately -- each entry is a claim that the action is verified,
        # made here because the check itself cannot reach the Marketplace.
        marketplace_verified = ("oven-sh/setup-bun",)

        violations = String[]
        workflows = sort(readdir(joinpath(REPO_ROOT, ".github", "workflows")))
        @test !isempty(workflows)
        for file in workflows
            (endswith(file, ".yml") || endswith(file, ".yaml")) || continue
            workflow = YAML.load_file(joinpath(REPO_ROOT, ".github", "workflows", file))
            for (_, job) in get(workflow, "jobs", Dict())
                for step in get(job, "steps", [])
                    uses = get(step, "uses", nothing)
                    uses === nothing && continue
                    startswith(uses, "./") && continue
                    parts = split(uses, '@'; limit = 2)
                    if length(parts) != 2
                        push!(violations, "$file: $uses has no @ref")
                        continue
                    end
                    src, ref = parts
                    if !occursin(r"^[0-9a-f]{40}$", ref)
                        push!(violations, "$file: $uses is not pinned to a full-length commit SHA")
                    end
                    owner = first(split(src, '/'))
                    if !(owner in github_owned_owners || src in marketplace_verified)
                        push!(violations, "$file: $uses is neither github-owned, hyperpolymath-owned, nor on the verified list")
                    end
                end
            end
        end
        @test isempty(violations)
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
        # MEASURED 2026-09-21 on PR #39: a job named exactly `Julia tests` with a 1x1
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
