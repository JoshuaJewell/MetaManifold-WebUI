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
const INSTALL_JL  = joinpath(REPO_ROOT, "scripts", "install.jl")
const CI_PATH     = joinpath(REPO_ROOT, ".github", "workflows", "ci.yml")
const WORKFLOWS   = joinpath(REPO_ROOT, ".github", "workflows")
const MANIFEST    = joinpath(REPO_ROOT, "Manifest.toml")
const RENV_LOCK   = joinpath(REPO_ROOT, "renv.lock")

# The tools the pipeline shells out to. A tool absent from the pin file is a tool
# nobody has decided the version of.
const EXPECTED_TOOLS = ["cutadapt", "fastqc", "multiqc", "vsearch", "cd_hit_est", "swarm",
                        "mafft", "trimal", "iqtree", "raxml", "gappa"]

# The platforms install.jl can resolve an archive for, as "<OS_TYPE>-<ARCH_STR>".
const EXPECTED_PLATFORMS = ["linux-x86_64", "linux-aarch64", "macos-x86_64", "macos-aarch64"]

is_sha256(s) = s isa AbstractString && occursin(r"^[0-9a-f]{64}$", s)

# This fork's account allows only GitHub's own actions and its own, each pinned to a
# full commit SHA.
allowed_action(uses) =
    occursin(r"^(actions|hyperpolymath)/[A-Za-z0-9_.-]+(/[A-Za-z0-9_./-]+)?@[0-9a-f]{40}$", uses)


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

        # Julia is the one runtime CI cannot read from the pin file, because nothing can
        # be read before Julia exists. Its archive URL and checksum are therefore
        # duplicated in the job's env, and this is the test that makes the duplication
        # safe. The pinned URL names the version, which must equal Manifest.toml's
        # (tested above).
        #
        # Read from the job's env rather than from a matrix: the matrix was removed
        # because GitHub appends a matrix combination to the posted check name (see "a
        # required check name is a stable identifier" below).
        archive = pins["runtimes"]["julia"]["archives"]["linux-x86_64"]
        @test job["env"]["JULIA_URL"] == archive["url"]
        @test job["env"]["JULIA_SHA256"] == archive["sha256"]

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
        for ref in (raw"$VSEARCH_URL", raw"$SWARM_URL", raw"$FASTQC_URL",
                    raw"$JULIA_URL", raw"$BUN_URL")
            @test occursin("fetch_pinned \"$ref\"", runs)
        end

        # Julia and Bun, which actions used to install, are refused unless they hash
        # to their pins.
        for (sha, file) in ((raw"$JULIA_SHA256", "julia.tar.gz"), (raw"$BUN_SHA256", "bun.zip"))
            @test occursin("echo \"$sha  $file\" | sha256sum -c -", runs)
        end
    end

    @testset "every action is one this account allows" begin
        # Anything else ends the whole workflow in startup_failure before any job runs,
        # so the run names no failing step. That is how ci.yml failed on every run on
        # this fork, from at least 2026-10-02 until it stopped using julia-actions,
        # oven-sh/setup-bun and codecov/codecov-action.
        files = filter(f -> endswith(f, ".yml") || endswith(f, ".yaml"), readdir(WORKFLOWS))
        @test "ci.yml" in files
        @testset "$file" for file in files
            for (_, job) in YAML.load_file(joinpath(WORKFLOWS, file))["jobs"],
                step in get(job, "steps", [])
                uses = get(step, "uses", nothing)
                uses === nothing || @test allowed_action(uses)
            end
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

@testset "the root keeps its shape" begin
    # Adapted from the RSR root-cleanup doctrine (rsr-template-repo,
    # docs/proposals/root-cleanup.adoc): the repository root holds authority
    # files, the entry points a newcomer runs, tool-mandated build metadata, and
    # one directory per concern -- nothing else. The check is BIDIRECTIONAL,
    # because a one-directional allowlist only ever ratchets open: a tracked
    # root entry that is not listed here is drift, and a listed entry that has
    # gone missing is drift too. The tracked tree is what is checked, so local
    # untracked build output (bin/, web/) does not fail a developer's run. To
    # add a root entry, add it here with its reason in the same commit.
    root_entries = [
        # Authority files. GitHub reads CITATION.cff from the default-branch
        # root only, so it cannot move.
        "README.md", "LICENSE", "CITATION.cff",
        # The two entry points a newcomer runs (the launcher stays at the
        # root); the installer itself lives under scripts/.
        "install.sh", "start.sh",
        # Build metadata rooted by its tooling: Julia Pkg, and renv, which
        # also scans R/ for dependency discovery.
        "Project.toml", "Manifest.toml", "renv.lock", "renv", ".Rprofile", "R",
        # One directory per concern. Health files (CONTRIBUTING.md) and CI
        # configs (codecov.yml) live under .github/.
        ".github", "bench", "config", "data", "docs", "frontend",
        "proofs", "scripts", "src", "test",
        # Agent instructions (the agents.md convention): coding agents load
        # this at startup, and it only works from the root, which is why it
        # is listed here rather than living under .github/ with the other
        # health files.
        "AGENTS.md",
        # VCS housekeeping.
        ".gitignore",
    ]
    tracked = Set(split(readchomp(`git -C $(REPO_ROOT) ls-tree --name-only HEAD`), "\n"))
    for entry in root_entries
        @test entry in tracked
    end
    extras = sort(collect(setdiff(tracked, Set(root_entries))))
    if !isempty(extras)
        println("root entries outside the allowlist: ", join(extras, ", "))
    end
    @test isempty(extras)
end

@testset "dependabot never batches a major" begin
    # #124 grouped every frontend update, and with nothing requiring the
    # checks it merged five framework majors (React 18->19, vite 6->8,
    # TypeScript 5->7, react-router-dom 6->7, plotly 2->4) in one PR, with
    # tsc already red until #140 repaired it. A major is a migration, not an
    # update, so it must open its own PR: a Dependabot group only claims an
    # update when its patterns match AND its update-types admit it -- a
    # dependency matching no group gets its own PR. Every group must
    # therefore NAME its update-types, and never "major". The frontend must
    # also stay on the "bun" ecosystem: the 2026-09 security-update runs
    # proved this account's npm_and_yarn updater refuses a Bun-managed
    # project (see the header of .github/dependabot.yml).
    dependabot = YAML.load_file(joinpath(REPO_ROOT, ".github", "dependabot.yml"))
    @test dependabot["version"] == 2
    ecosystems = Dict{String,Any}()
    for u in dependabot["updates"]
        @test !haskey(ecosystems, u["package-ecosystem"])
        ecosystems[u["package-ecosystem"]] = u
    end
    @test sort(collect(keys(ecosystems))) == ["bun", "github-actions"]
    @test ecosystems["bun"]["directory"] == "/frontend"
    for u in dependabot["updates"]
        for (_, group) in get(u, "groups", Dict{String,Any}())
            update_types = get(group, "update-types", nothing)
            @test update_types !== nothing
            @test !("major" in update_types)
        end
    end
end
