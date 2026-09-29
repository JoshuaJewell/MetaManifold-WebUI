#!/usr/bin/env julia
#
# Dependency installer for MetaManifold
#
# Installs Julia deps, checks/downloads external CLI tools, and installs
# required R packages. Writes resolved tool paths to config/tools.yml.
#
# Usage via install.sh, or directly by:
#   julia --project=. install.jl [--update] [--modify] [--sysimage]
#
# Options:
#   --update    Re-assert the pinned tool versions and refresh managed binaries in bin/
#   --modify    Revisit configured tool paths instead of silently reusing them
#   --sysimage  After installation, compile a sysimage of all Julia deps to
#               speed up subsequent startup. Output: MetaManifold.so (.dylib on macOS)
#               Use with: julia --sysimage MetaManifold.so --project=. ...
#
# Pinning: every external tool version, download URL, and archive checksum lives in
# config/defaults/tool_versions.yml, and nothing in this script tracks "latest".
# Two clean installs a year apart therefore obtain the same binaries. --update
# re-fetches the same pins, so it is idempotent; a version moves only when that
# file is edited. An archive whose SHA256 does not match its pin is refused.

using Pkg

# False when the test suite loads this file to exercise the pin lookup and the
# checksum refusal. Nothing may be installed, and no dependency resolved, on a load
# that is not an install.
const RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == abspath(@__FILE__)

if RUNNING_AS_SCRIPT
    @info "Installing Julia package dependencies..."
    Pkg.instantiate()
end

using YAML
using SHA
import Downloads

## Instantiate
const UPDATE_MODE   = "--update"   in ARGS
const MODIFY_MODE   = "--modify"   in ARGS
const SYSIMAGE_MODE = "--sysimage" in ARGS
const PROJECT_ROOT  = @__DIR__
const BIN_DIR       = joinpath(PROJECT_ROOT, "bin")
const CONFIG_DIR    = joinpath(PROJECT_ROOT, "config")
const TOOLS_CONFIG  = joinpath(CONFIG_DIR, "tools.yml")
const VERSIONS_FILE = joinpath(CONFIG_DIR, "defaults", "tool_versions.yml")

const OS_TYPE = Sys.islinux() ? "linux" :
                Sys.isapple() ? "macos" :
                error("Unsupported OS. Only Linux and macOS are supported.")

const ARCH_STR = Sys.ARCH == :x86_64  ? "x86_64"  :
                 Sys.ARCH == :aarch64 ? "aarch64" :
                 string(Sys.ARCH)

# Canonical binary name for tools whose config key differs from the binary name.
const BINARY_NAMES = Dict("cd_hit_est" => "cd-hit-est", "iqtree" => "iqtree3",
                          "raxml" => "raxmlHPC-PTHREADS-SSE3")
bin_name(key::String) = get(BINARY_NAMES, key, key)

mkpath(BIN_DIR)

## Install summary
#
# Every step that resolves, defers, or fails a dependency records one line here,
# and main() prints the collected block last. The rule this enforces: nothing the
# installer chose not to do may pass without a line the operator can read. A
# skipped R runtime, an unresolved tool, headers that need root - each is reported
# at install time.

@enum StepStatus STEP_OK STEP_SKIPPED STEP_FAILED STEP_ACTION

const SUMMARY = Tuple{String,StepStatus,String}[]

record!(name::AbstractString, status::StepStatus, detail::AbstractString = "") =
    push!(SUMMARY, (String(name), status, String(detail)))

function print_summary()
    isempty(SUMMARY) && return
    labelw = maximum(length(s[1]) for s in SUMMARY)
    println()
    println("=================== Install summary ===================")
    for (name, status, detail) in SUMMARY
        tag = status == STEP_OK      ? "OK"            :
              status == STEP_SKIPPED ? "SKIPPED"       :
              status == STEP_FAILED  ? "FAILED"        :
                                       "ACTION NEEDED"
        line = "  $(rpad(name, labelw))  $(rpad(tag, 13))"
        println(isempty(detail) ? rstrip(line) : "$line $detail")
    end

    pending = filter(s -> s[2] != STEP_OK, SUMMARY)
    if isempty(pending)
        println()
        println("  Everything the installer manages is in place. Start with:  bash start.sh")
    else
        println()
        println("  Still needs attention:")
        for (name, _, detail) in pending
            println(isempty(detail) ? "    - $name" : "    - $name -> $detail")
        end
    end
    println("======================================================")
end

## Version pins
# The pin file is the sole source of truth for what this installer fetches. Without
# it there is no defensible version to install, so its absence is fatal rather than
# an invitation to fall back on whatever upstream published most recently.
isfile(VERSIONS_FILE) || error(
    "Version pin file not found: $VERSIONS_FILE\n" *
    "It is committed to the repository; a checkout missing it is incomplete."
)
const PINS     = YAML.load_file(VERSIONS_FILE)
const PLATFORM = "$(OS_TYPE)-$(ARCH_STR)"

pinned_version(tool::String)::String = string(PINS["tools"][tool]["version"])

# The archive pinned for this platform, or the "any" entry for artefacts that are
# not platform-specific (FastQC is Java; cd-hit is compiled from source).
function pinned_archive(tool::String)::Dict
    archives = PINS["tools"][tool]["archives"]
    rec = get(archives, PLATFORM, get(archives, "any", nothing))
    rec === nothing && error(
        "No $tool archive is pinned for $PLATFORM in $VERSIONS_FILE.\n" *
        "Add one, with its SHA256, or install $tool yourself and give install.jl the path."
    )
    rec
end

## Config loading / saving
function load_tools_config()::Dict{String,Any}
    isfile(TOOLS_CONFIG) || return Dict{String,Any}()
    data = YAML.load_file(TOOLS_CONFIG)
    data isa Dict ? data : Dict{String,Any}()
end

function write_tools_config(config::Dict)
    mkpath(CONFIG_DIR)
    escaped = Dict{String,Any}()
    for key in sort(collect(keys(config)))
        val = config[key]
        path = val isa Dict ? get(val, "path", nothing) : val
        escaped[key] = Dict("path" => path)
    end

    yaml = YAML.write(escaped)
    open(TOOLS_CONFIG, "w") do io
        print(io, yaml)
    end
    @info "Config written to $TOOLS_CONFIG"
end

## Checking for tool
# Returns true if the binary at `path` is callable (local or remote SSH).
function check_tool(path::String)::Bool
    if occursin('@', path)
        # Remote SSH path: user@host:/path/to/binary
        colon_idx = findfirst(':', path)
        colon_idx === nothing && return false
        user_host = path[1:colon_idx-1]
        bin_path  = path[colon_idx+1:end]
        try
            run(pipeline(
                `ssh -o BatchMode=yes -o ConnectTimeout=5 $user_host test -x $bin_path`;
                stdout=devnull, stderr=devnull
            ))
            return true
        catch
            return false
        end
    else
        # Local path or bare name (PATH lookup)
        resolved = isfile(path) ? path : Sys.which(path)
        resolved === nothing && return false
        return Sys.isexecutable(resolved)
    end
end

"""Return the full path to `name` if it is in PATH, otherwise nothing."""
function find_in_path(name::String)
    Sys.which(name)
end

## Interactive prompts
function prompt_yn(question::String, default_yes::Bool = true)::Bool
    hint = default_yes ? "[Y/n]" : "[y/N]"
    print("  $question $hint: ")
    answer = strip(readline())
    isempty(answer) && return default_yes
    return lowercase(answer) in ("y", "yes")
end

function prompt_path(label::String)::Union{String,Nothing}
    print("  $label: ")
    p = strip(readline())
    isempty(p) ? nothing : p
end

## Download helpers
# Download `url` to `dest`.
function download_to(url::String, dest::String)
    @info "Downloading $(basename(url))..."
    Downloads.download(url, dest)
end

# Recursively check `dir` and return the first file named `filename`, or nothing.
function find_file_in_dir(dir::String, filename::String)::Union{String,Nothing}
    for (root, _dirs, files) in walkdir(dir)
        idx = findfirst(==(filename), files)
        idx !== nothing && return joinpath(root, files[idx])
    end
    nothing
end

# Reject the file unless it hashes to `expected`. A missing pin also aborts. The
# offending file is deleted so a failed run leaves no unverified artefact behind.
function verify_sha256(path::String, expected, source::String)
    expected === nothing && error(
        "No SHA256 is pinned for $source in $VERSIONS_FILE.\n" *
        "Record the checksum before installing; an unverifiable download is refused."
    )
    actual = bytes2hex(open(sha256, path))
    want   = lowercase(strip(string(expected)))
    if actual != want
        rm(path; force=true)
        error(
            "Checksum mismatch for $(basename(source)).\n" *
            "  expected: $want\n" *
            "  actual:   $actual\n" *
            "The download does not match its pin in $VERSIONS_FILE, so it is not the\n" *
            "artefact this project was tested against. Refusing to install it."
        )
    end
    @info "Verified SHA256 of $(basename(source))."
end

function download_verified(url::String, dest::String, expected)
    download_to(url, dest)
    verify_sha256(dest, expected, url)
    dest
end

## Tool download functions
# Unpack `tarball` in a scratch directory and move `binary` out of it into bin/.
# Unpacking away from bin/ is what keeps the result deterministic: a source tree
# left by an earlier version can otherwise be picked up in place of this one.
function extract_binary(tarball::String, binary::String)::String
    workdir = mktempdir()
    try
        if endswith(tarball, ".zip")
            run(`unzip -q -o $tarball -d $workdir`)
        else
            run(`tar -xzf $tarball -C $workdir --warning=no-unknown-keyword`)
        end
        found = find_file_in_dir(workdir, binary)
        found === nothing && error("$binary not found in $(basename(tarball)) after extraction.")
        dest = joinpath(BIN_DIR, binary)
        mv(found, dest; force=true)
        chmod(dest, 0o755)
        return dest
    finally
        rm(workdir; recursive=true, force=true)
    end
end

function download_pinned_binary(tool::String, binary::String)::String
    rec     = pinned_archive(tool)
    version = pinned_version(tool)
    @info "Installing $tool $version (pinned)..."

    ext     = endswith(string(rec["url"]), ".zip") ? ".zip" : ".tar.gz"
    tarball = joinpath(BIN_DIR, "$(binary)_download$ext")
    try
        download_verified(rec["url"], tarball, get(rec, "sha256", nothing))
        return extract_binary(tarball, binary)
    finally
        rm(tarball; force=true)
    end
end

download_vsearch()::String = download_pinned_binary("vsearch", "vsearch")
download_trimal()::String  = download_pinned_binary("trimal", "trimal")
download_iqtree()::String  = download_pinned_binary("iqtree", "iqtree3")

# Unpack a pinned source tarball in a scratch directory and hand its top-level
# directory to `build`, which returns the built binary's path.
function build_pinned_source(build::Function, tool::String)::String
    rec     = pinned_archive(tool)
    version = pinned_version(tool)
    Sys.which("make") === nothing && error(
        "'make' is not available to build $tool $version from source.\n" *
        "Install $tool yourself and enter its path when prompted.")
    _c_compiler() === nothing && error(
        "No C compiler is available to build $tool $version from source.\n" *
        "Install build tools (e.g. build-essential), or install $tool yourself.")
    tarball = joinpath(BIN_DIR, "$(tool)_source.tar.gz")
    workdir = mktempdir()
    try
        @info "Building $tool $version from its pinned source..."
        download_verified(string(rec["url"]), tarball, get(rec, "sha256", nothing))
        run(`tar -xzf $tarball -C $workdir --warning=no-unknown-keyword`)
        dirs = filter(isdir, readdir(workdir; join=true))
        length(dirs) == 1 || error("Unexpected layout in the $tool source tarball.")
        return build(only(dirs))
    finally
        rm(tarball; force=true)
        rm(workdir; recursive=true, force=true)
    end
end

# A distribution package is whatever that distribution ships, which need not be
# the pinned version; the version actually obtained is recorded at preflight.
function try_package(pkg::String, binaries::Vector{String}; brew_pkg=pkg)
    cmd = package_install_cmd([pkg]; brew_pkg)
    cmd === nothing && return nothing
    @info "Trying package manager install for $pkg..."
    try
        run(cmd)
    catch
        @warn "Package manager install of $pkg failed, falling back to a source build."
        return nothing
    end
    for b in binaries
        found = Sys.which(b)
        found !== nothing && return found
    end
    nothing
end

function download_mafft()::String
    found = try_package("mafft", ["mafft"])
    found !== nothing && return found
    prefix = joinpath(BIN_DIR, "mafft")
    build_pinned_source("mafft") do src
        run(Cmd(`make -j$(Sys.CPU_THREADS) PREFIX=$prefix install`; dir=joinpath(src, "core")))
        bin = joinpath(prefix, "bin", "mafft")
        isfile(bin) || error("mafft not found at $bin after building.")
        bin
    end
end

# The PTHREADS build is what -T needs; SSE3 exists only on x86.
function download_raxml()::String
    found = try_package("raxml", ["raxmlHPC-PTHREADS-SSE3", "raxmlHPC-PTHREADS"])
    found !== nothing && return found
    sse = ARCH_STR == "x86_64"
    binary   = sse ? "raxmlHPC-PTHREADS-SSE3" : "raxmlHPC-PTHREADS"
    makefile = "Makefile." * (sse ? "SSE3." : "") * "PTHREADS." * (OS_TYPE == "macos" ? "mac" : "gcc")
    build_pinned_source("raxml") do src
        run(Cmd(`make -f $makefile`; dir=src))
        bin = joinpath(src, binary)
        isfile(bin) || error("$binary not found after building RAxML.")
        dest = joinpath(BIN_DIR, binary)
        cp(bin, dest; force=true)
        chmod(dest, 0o755)
        dest
    end
end

# gappa's build downloads genesis, CLI11 and sparsepp at the commits its
# CMakeLists names, so it needs cmake, a C++17 compiler and network access.
function download_gappa()::String
    Sys.which("cmake") === nothing && error(
        "gappa is built from source and needs cmake. Install cmake and a C++ compiler, " *
        "or build gappa yourself and enter the path to bin/gappa.")
    build_pinned_source("gappa") do src
        run(Cmd(`make`; dir=src))
        bin = joinpath(src, "bin", "gappa")
        isfile(bin) || error("gappa not found at $bin after building.")
        dest = joinpath(BIN_DIR, "gappa")
        cp(bin, dest; force=true)
        chmod(dest, 0o755)
        dest
    end
end

function download_fastqc()::String
    rec     = pinned_archive("fastqc")
    version = pinned_version("fastqc")
    url     = string(rec["url"])

    zipfile = joinpath(BIN_DIR, "fastqc.zip")
    try
        download_to(url, zipfile)
    catch e
        rm(zipfile; force=true)
        error(
            "Failed to download FastQC v$version from $url\n" *
            "Babraham publish no release API, so this URL is pinned by hand in\n" *
            "$VERSIONS_FILE and may have moved. Download it manually from\n" *
            "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ and place the\n" *
            "fastqc binary in bin/.\n" *
            "Original error: $e"
        )
    end
    verify_sha256(zipfile, get(rec, "sha256", nothing), url)

    run(`unzip -q -o $zipfile -d $BIN_DIR`)
    rm(zipfile)

    bin = joinpath(BIN_DIR, "FastQC", "fastqc")
    isfile(bin) || error("fastqc not found after extraction. Expected at $bin")
    chmod(bin, 0o755)
    bin
end

## FastQC's Java runtime
# FastQC ships as a Perl wrapper around a Java application, so a JRE is a hard
# runtime dependency that downloading FastQC cannot itself satisfy. Left
# unchecked the failure is near-silent: the wrapper writes "Can't exec java" to
# stderr, prints nothing to stdout and still exits 0, so the shortfall surfaces
# much later as an empty version string from the provenance probe, or as a
# fastqc stage that fails without ever naming Java. Naming it here, with the
# command that fixes it, keeps the diagnosis at install time.

# The install command to advise the user to run. It is returned even without root,
# since the user runs it themselves.
function java_install_hint()::String
    OS_TYPE == "macos" && return "brew install openjdk"
    prefix = has_root() ? "" : "sudo "
    Sys.which("apt")    !== nothing && return "$(prefix)apt install -y default-jre-headless"
    Sys.which("dnf")    !== nothing && return "$(prefix)dnf install -y java-17-openjdk-headless"
    Sys.which("pacman") !== nothing && return "$(prefix)pacman -S --needed --noconfirm jre-openjdk-headless"
    Sys.which("zypper") !== nothing && return "$(prefix)zypper install -y java-17-openjdk-headless"
    "install a Java runtime (JRE 11 or newer) using your system package manager"
end

# Verify FastQC can actually run, by the same measure the provenance probe uses:
# a version banner on stdout. The wrapper's exit status is not that measure --
# it is 0 whether or not Java was found.
function check_fastqc_runtime(fastqc::AbstractString)
    out, err = IOBuffer(), IOBuffer()
    try
        run(pipeline(`$fastqc --version`; stdout=out, stderr=err))
    catch
        # A non-zero exit is just another way to fail; the banner check below decides.
    end
    banner = strip(String(take!(out)))
    complaint = strip(String(take!(err)))

    if occursin(r"FastQC\s+v"i, banner)
        record!("Java runtime", STEP_OK, isempty(banner) ? "present" : banner)
        return true
    end

    if Sys.which("java") === nothing
        hint = java_install_hint()
        @warn "FastQC is installed but cannot run: no Java runtime on PATH.\n" *
              "FastQC is a Java application, so the fastqc stage will fail until a JRE\n" *
              "is installed. Fix it with:\n\n    $hint\n"
        record!("Java runtime", STEP_ACTION,
                "missing - FastQC cannot run without it. Install with:  $hint")
        return false
    end

    # Java is present but FastQC still would not report a version: surface what it said.
    detail = isempty(complaint) ? "no version banner from `$fastqc --version`" :
                                  first(split(complaint, '\n'))
    @warn "FastQC is installed and Java is on PATH, but FastQC did not report a version.\n" *
          "The fastqc stage may fail. It said: $detail"
    record!("Java runtime", STEP_ACTION, "Java found, but FastQC failed: $detail")
    false
end

function download_cdhit()::String
    version = pinned_version("cd_hit_est")

    pkg_cmd = package_install_cmd(["cd-hit"]; brew_pkg="cd-hit")
    if pkg_cmd !== nothing
        # A distribution package is whatever that distribution ships, which need not
        # be the pinned version. It is preferred anyway, because building cd-hit from
        # source needs a C++ toolchain that many machines lack; the version actually
        # obtained is recorded at preflight, where a discrepancy is visible.
        @info "Trying package manager install for cd-hit (expected v$version)..."
        try
            run(pkg_cmd)
            found = Sys.which("cd-hit-est")
            found !== nothing && return found
        catch
            @warn "Package manager install failed, falling back to source build."
        end
    end

    # Fallback: build from the pinned source tarball (upstream publish no precompiled
    # Linux binary).
    Sys.which("make") === nothing && error(
        "cd-hit could not be installed via package manager and 'make' is not available to build from source.\n" *
        "Install cd-hit manually and enter the path to cd-hit-est when prompted."
    )

    rec     = pinned_archive("cd_hit_est")
    tarball = joinpath(BIN_DIR, "cdhit_download.tar.gz")
    workdir = mktempdir()
    try
        @info "Downloading cd-hit v$version source and building from source..."
        download_verified(string(rec["url"]), tarball, get(rec, "sha256", nothing))
        run(`tar -xzf $tarball -C $workdir --warning=no-unknown-keyword`)

        src_dir = nothing
        for entry in readdir(workdir; join=true)
            isdir(entry) && startswith(basename(entry), "cd-hit") && (src_dir = entry; break)
        end
        src_dir === nothing && error("cd-hit source directory not found after extraction.")

        run(Cmd(`make -j$(Sys.CPU_THREADS)`; dir=src_dir))

        bin = joinpath(src_dir, "cd-hit-est")
        isfile(bin) || error("cd-hit-est binary not found after building. Check that a C++ compiler is installed.")

        dest = joinpath(BIN_DIR, "cd-hit-est")
        cp(bin, dest; force=true)
        chmod(dest, 0o755)
        return dest
    finally
        rm(tarball; force=true)
        rm(workdir; recursive=true, force=true)
    end
end

"""
Add the user's local bin directory to ENV["PATH"] so that Sys.which() and
subsequent Cmd calls can find freshly-installed scripts without restarting.
"""
function ensure_local_bin_on_path()
    local_bin = joinpath(homedir(), ".local", "bin")
    paths = split(get(ENV, "PATH", ""), ':')
    if local_bin ∉ paths
        ENV["PATH"] = local_bin * ":" * ENV["PATH"]
    end
end

function ensure_pipx()
    Sys.which("pipx") !== nothing && return  # already present

    @info "pipx not found - attempting to install it..."

    pkg_cmd = package_install_cmd(["pipx"]; brew_pkg="pipx", pacman_pkg="python-pipx", zypper_pkg="python3-pipx")
    if pkg_cmd !== nothing
        try
            run(pkg_cmd)
            run(`pipx ensurepath`)
            ensure_local_bin_on_path()
            Sys.which("pipx") !== nothing && return
        catch
            @warn "Package manager install of pipx failed."
        end
    end

    # Fallback: bootstrap pipx via pip/python3
    pip_cmd = nothing
    for candidate in (`pip3`, `pip`, `python3 -m pip`)
        try
            run(pipeline(`$candidate --version`; stdout=devnull, stderr=devnull))
            pip_cmd = candidate
            break
        catch
        end
    end

    if pip_cmd !== nothing
        try
            run(`$pip_cmd install --user pipx`)
            run(`python3 -m pipx ensurepath`)
            # Update PATH in the running process so Sys.which finds pipx
            ensure_local_bin_on_path()
            Sys.which("pipx") !== nothing && return
        catch
        end
    end

    @warn "Could not install pipx automatically. Python tools (cutadapt, multiqc) " *
          "may need to be installed manually."
end

function pipx_has_tool(name::String)::Bool
    Sys.which("pipx") === nothing && return false
    try
        output = read(`pipx list --short`, String)
        # Each line reads "<package> <version>", so only the first field is the name.
        for line in split(output, '\n')
            fields = split(strip(line))
            !isempty(fields) && fields[1] == name && return true
        end
        return false
    catch
        return false
    end
end

function install_python_tool(name::String)::String
    version = pinned_version(name)
    spec    = "$(name)==$(version)"

    # pipx (recommended on PEP 668 / Debian-managed systems)
    if Sys.which("pipx") !== nothing
        # There is no upgrade path, by design: `pipx upgrade` would walk the tool off
        # its pin. Where the tool is already present, --force reinstalls it at the
        # pinned version, which is also what makes --update idempotent rather than
        # a slow drift towards whatever PyPI published last.
        args = pipx_has_tool(name) ?
            ["pipx", "install", "--force", spec] :
            ["pipx", "install", spec]

        @info "Installing $spec via pipx..."
        run(Cmd(args))
        ensure_local_bin_on_path()
        found = find_in_path(name)
        found !== nothing && return found
        # pipx installs to ~/.local/bin by default
        local_bin = joinpath(homedir(), ".local", "bin", name)
        isfile(local_bin) && return local_bin
        @warn "Installed $name via pipx but could not locate the binary. Ensure ~/.local/bin is in PATH."
        return name
    end

    # Fallback to pip --user, then --break-system-packages if blocked
    pip_cmd = nothing
    for candidate in (`pip3`, `pip`, `python3 -m pip`)
        try
            run(pipeline(`$candidate --version`; stdout=devnull, stderr=devnull))
            pip_cmd = candidate
            break
        catch
        end
    end
    pip_cmd === nothing && error(
        "Neither pipx nor pip found. Install pipx (recommended) or Python 3 with pip."
    )

    @info "Installing $spec via pip..."
    success = try
        run(`$pip_cmd install --user $spec`)
        true
    catch
        false
    end

    if !success
        # Fallback 2 to PEP 668: externally-managed environment, try --break-system-packages
        @warn "pip --user blocked by system policy. Retrying with --break-system-packages..."
        run(`$pip_cmd install --user --break-system-packages $spec`)
    end

    ensure_local_bin_on_path()

    # Locate the installed binary
    found = find_in_path(name)
    found !== nothing && return found

    local_bin = joinpath(homedir(), ".local", "bin", name)
    isfile(local_bin) && return local_bin

    @warn "Installed $name but could not locate the binary. Ensure ~/.local/bin is in PATH."
    name
end

function download_swarm()::String
    # An already-present bin/swarm is honoured, but not under --update, whose whole
    # purpose is to re-assert the pin. Trusting it there would let a binary of
    # unknown provenance, predating the pins, survive every update indefinitely.
    bundled = joinpath(BIN_DIR, "swarm")
    if !UPDATE_MODE && isfile(bundled) && Sys.isexecutable(bundled)
        @info "Using bundled swarm binary at $bundled"
        return bundled
    end

    download_pinned_binary("swarm", "swarm")
end

function install_r_sysdeps()
    # System libraries required to compile Bioconductor / tidyverse packages from source.
    # Package names differ across distro families; each list maps to the same underlying
    # libraries (bzip2, xz, zlib, curl, openssl, libxml2, freetype, libpng, libjpeg,
    # libtiff, fontconfig, harfbuzz, fribidi, hdf5).

    apt_deps = [
        "pkg-config",
        "libbz2-dev", "liblzma-dev", "zlib1g-dev",
        "libcurl4-openssl-dev", "libssl-dev",
        "libxml2-dev",
        "libfreetype6-dev", "libpng-dev",
        "libjpeg-dev", "libtiff5-dev",
        "libfontconfig1-dev",
        "libharfbuzz-dev", "libfribidi-dev",
        "libhdf5-dev",
    ]

    dnf_deps = [
        "pkgconf-pkg-config",
        "bzip2-devel", "xz-devel", "zlib-devel",
        "libcurl-devel", "openssl-devel",
        "libxml2-devel",
        "freetype-devel", "libpng-devel",
        "libjpeg-turbo-devel", "libtiff-devel",
        "fontconfig-devel",
        "harfbuzz-devel", "fribidi-devel",
        "hdf5-devel",
    ]

    pacman_deps = [
        "pkgconf",
        "bzip2", "xz", "zlib",
        "curl", "openssl",
        "libxml2",
        "freetype2", "libpng",
        "libjpeg-turbo", "libtiff",
        "fontconfig",
        "harfbuzz", "fribidi",
        "hdf5",
    ]

    zypper_deps = [
        "pkg-config",
        "libbz2-devel", "xz-devel", "zlib-devel",
        "libcurl-devel", "libopenssl-devel",
        "libxml2-devel",
        "freetype2-devel", "libpng16-devel",
        "libjpeg8-devel", "libtiff-devel",
        "fontconfig-devel",
        "harfbuzz-devel", "fribidi-devel",
        "hdf5-devel",
    ]

    OS_TYPE != "linux" && return

    if Sys.which("apt") !== nothing
        _install_sysdeps_with(package_install_cmd(apt_deps; linux_manager=:apt), apt_deps, "apt")
    elseif Sys.which("dnf") !== nothing
        _install_sysdeps_with(package_install_cmd(dnf_deps; linux_manager=:dnf), dnf_deps, "dnf")
    elseif Sys.which("pacman") !== nothing
        _install_sysdeps_with(package_install_cmd(pacman_deps; linux_manager=:pacman), pacman_deps, "pacman")
    elseif Sys.which("zypper") !== nothing
        _install_sysdeps_with(package_install_cmd(zypper_deps; linux_manager=:zypper), zypper_deps, "zypper")
    else
        @warn "Could not detect a supported package manager (apt, dnf, pacman, zypper). " *
              "Some R packages may fail to compile. Install the development headers for: " *
              "bzip2, xz, zlib, curl, openssl, libxml2, freetype, libpng, libjpeg, " *
              "libtiff, fontconfig, harfbuzz, fribidi, hdf5"
    end
end

function _install_sysdeps_with(cmd::Union{Cmd,Nothing}, deps::Vector{String}, label::String)
    if cmd === nothing
        @info "Skipping automated $label install for R system dependencies because this session does not have root or passwordless sudo. " *
              "If compilation fails, install these packages manually:\n  " * join(deps, " ")
        return
    end
    @info "Installing R system library dependencies via $label..."
    try
        run(cmd)
    catch
        @warn "$label install of R system deps failed. " *
              "Some R packages may not compile; install these packages manually:\n  " * join(deps, " ")
    end
end

# The R side is renv-managed: renv.lock is the single source of truth and the
# runtime activates it through .Rprofile. The installer therefore reproduces the
# library straight from the lockfile with renv::restore(), which writes only to
# the project-local renv/library and never needs root. renv itself is
# self-bootstrapping from renv/activate.R, so a machine with a bare R can still
# run this.

# Build-time libraries the Bioconductor/tidyverse stack compiles against; a
# missing one is the usual reason renv::restore() dies partway (e.g. Rhtslib
# needs curl/curl.h). Each row: (probe, apt, dnf, pacman, zypper).
#
# `probe` is "pc:<module>" for a pkg-config query, or "h:<header>" for a compiler
# `#include` test. pkg-config is authoritative where a .pc file exists; the
# header test covers libraries that ship none (libbz2-dev on Debian has no .pc),
# and libraries whose headers sit on a non-default include path (libxml2,
# freetype2, harfbuzz, fribidi) are left as pc: since a bare #include would
# wrongly report them absent.
const R_SYSDEP_TABLE = [
    #  probe             apt                     dnf                pacman        zypper
    ("pc:libcurl",     "libcurl4-openssl-dev", "libcurl-devel",    "curl",       "libcurl-devel"),
    ("pc:openssl",     "libssl-dev",           "openssl-devel",    "openssl",    "libopenssl-devel"),
    ("pc:libxml-2.0",  "libxml2-dev",          "libxml2-devel",    "libxml2",    "libxml2-devel"),
    ("h:zlib.h",       "zlib1g-dev",           "zlib-devel",       "zlib",       "zlib-devel"),
    ("h:bzlib.h",      "libbz2-dev",           "bzip2-devel",      "bzip2",      "libbz2-devel"),
    ("h:lzma.h",       "liblzma-dev",          "xz-devel",         "xz",         "xz-devel"),
    ("pc:freetype2",   "libfreetype6-dev",     "freetype-devel",   "freetype2",  "freetype2-devel"),
    ("h:png.h",        "libpng-dev",           "libpng-devel",     "libpng",     "libpng16-devel"),
    ("h:tiff.h",       "libtiff5-dev",         "libtiff-devel",    "libtiff",    "libtiff-devel"),
    ("pc:fontconfig",  "libfontconfig1-dev",   "fontconfig-devel", "fontconfig", "fontconfig-devel"),
    ("pc:harfbuzz",    "libharfbuzz-dev",      "harfbuzz-devel",   "harfbuzz",   "harfbuzz-devel"),
    ("pc:fribidi",     "libfribidi-dev",       "fribidi-devel",    "fribidi",    "fribidi-devel"),
]

const _CC_CACHE = Ref{Union{String,Nothing,Missing}}(missing)
function _c_compiler()::Union{String,Nothing}
    if _CC_CACHE[] === missing
        _CC_CACHE[] = nothing
        for c in ("cc", "gcc", "clang")
            p = Sys.which(c)
            if p !== nothing
                _CC_CACHE[] = p
                break
            end
        end
    end
    _CC_CACHE[]
end

# true / false, or nothing when the probe itself cannot run (no pkg-config, no
# compiler) - the caller reports "nothing" as needs-attention rather than
# assuming success.
function _r_dep_present(probe::AbstractString)::Union{Bool,Nothing}
    kind, arg = split(probe, ':'; limit = 2)
    if kind == "pc"
        pc = Sys.which("pkg-config")
        pc === nothing && return nothing
        return try; success(`$pc --exists $arg`); catch; false; end
    else
        cc = _c_compiler()
        cc === nothing && return nothing
        return mktemp() do path, io
            write(io, "#include <$arg>\n")
            close(io)
            try
                run(pipeline(`$cc -xc -fsyntax-only $path`; stdout = devnull, stderr = devnull))
                true
            catch
                false
            end
        end
    end
end

"""Rows of R_SYSDEP_TABLE whose library is missing or unverifiable. Linux only."""
function probe_missing_r_headers()::Vector{NTuple{5,String}}
    OS_TYPE == "linux" || return NTuple{5,String}[]
    missing = NTuple{5,String}[]
    if Sys.which("pkg-config") === nothing
        push!(missing, ("pc:pkg-config", "pkg-config", "pkgconf-pkg-config", "pkgconf", "pkg-config"))
    end
    if _c_compiler() === nothing
        push!(missing, ("cc", "build-essential", "gcc", "base-devel", "gcc"))
    end
    for row in R_SYSDEP_TABLE
        _r_dep_present(row[1]) == true || push!(missing, row)
    end
    missing
end

"""One copy-pasteable command to install the missing dev packages, or "" if none."""
function r_sysdep_hint(missing::Vector{NTuple{5,String}})::String
    isempty(missing) && return ""
    col, pre = Sys.which("apt")     !== nothing ? (2, "sudo apt install -y")     :
               Sys.which("dnf")     !== nothing ? (3, "sudo dnf install -y")     :
               Sys.which("pacman")  !== nothing ? (4, "sudo pacman -S --needed") :
               Sys.which("zypper")  !== nothing ? (5, "sudo zypper install -y")  :
                                                  (2, "install dev headers:")
    join([pre; unique(row[col] for row in missing)], " ")
end

# Reproduce renv/library from renv.lock. Returns (:ok, "") or (:incomplete, "pkg,pkg").
function setup_r_packages(; rebuild::Bool=false)::Tuple{Symbol,String}
    restore_call = rebuild ? "renv::restore(prompt = FALSE, rebuild = TRUE)" :
                             "renv::restore(prompt = FALSE)"
    try
        run(Cmd(`Rscript -e $restore_call`; dir = PROJECT_ROOT))   # streams renv's progress
    catch e
        @warn "renv::restore() exited non-zero; verifying what landed anyway: $e"
    end

    # Verdict is whether the packages the pipeline loads are actually usable, not
    # renv's exit code: restore can fail on one leaf package and still leave a
    # working library, or "succeed" against a stale cache.
    probe = raw"""
        pkgs <- c("dada2", "vegan", "dplyr", "tibble")
        cat(paste(pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)],
                  collapse = ","))
    """
    bad = try
        strip(read(Cmd(`Rscript -e $probe`; dir = PROJECT_ROOT), String))
    catch
        "dada2,vegan,dplyr,tibble"
    end

    isempty(bad) ? (:ok, "") : (:incomplete, String(bad))
end

## Per-tool resolution

# Resolves a single tool interactively. Returns the resolved path string, or nothing
# if the user chose to skip. When `install_fn` is provided, offers an auto-install
# option as the first choice.
function resolve_tool(
    key::String,
    display_name::String,
    existing_config::Dict,
    install_fn::Union{Function,Nothing} = nothing
)::Union{String,Nothing}
    bin = bin_name(key)
    heading_printed = false
    function show_heading()
        if !heading_printed
            println()
            println("  --- $display_name -------------------------------------------------")
            heading_printed = true
        end
    end

    # Check existing config (skip in update mode for managed installs)
    existing = get(existing_config, key, nothing)
    existing_path = existing isa Dict ? get(existing, "path", nothing) : nothing

    if existing_path !== nothing && !UPDATE_MODE
        if check_tool(string(existing_path))
            if MODIFY_MODE
                show_heading()
                println("  Configured path: $existing_path")
                prompt_yn("  Use this?") && return string(existing_path)
            else
                return string(existing_path)
            end
        else
            show_heading()
            println("  Configured path: $existing_path")
            println("  Warning: configured path does not appear to be callable.")
        end
    end

    # Check PATH
    path_result = find_in_path(bin)
    if path_result !== nothing
        show_heading()
        println("  Found in PATH:   $path_result")
        prompt_yn("  Use this?") && return path_result
    end

    # Build option list
    show_heading()
    options = String[]
    if install_fn !== nothing
        push!(options, "Install/download automatically to bin/")
    end
    push!(options, "Enter a path manually  (/path/to/$bin)")
    push!(options, "Skip  (configure later in config/tools.yml)")

    for (i, opt) in enumerate(options)
        println("  $i) $opt")
    end

    print("  Choice [1]: ")
    raw = strip(readline())
    choice = isempty(raw) ? 1 : something(tryparse(Int, raw), 1)

    if install_fn !== nothing
        if choice == 1
            try
                return install_fn()
            catch e
                @error "Auto-install failed: $e"
                println()
                println("  What would you like to do?")
                println("  1) Enter a path manually  (/path/to/$bin)")
                println("  2) Skip  (configure later in config/tools.yml)")
                print("  Choice [1]: ")
                raw2 = strip(readline())
                choice = isempty(raw2) ? 1 : something(tryparse(Int, raw2), 1)
                choice == 1 || return nothing
                return prompt_path("  Enter path")
            end
        end
        if choice == 2
            return prompt_path("  Enter path")
        end
        return nothing  # Skip
    else
        if choice == 1
            return prompt_path("  Enter path")
        end
        return nothing  # Skip
    end
end

## Frontend
# bun builds the frontend into web/dist, which the server serves. Only the build
# needs it; pipeline runs never do.
const FRONTEND_DIR = joinpath(PROJECT_ROOT, "frontend")
const DIST_INDEX   = joinpath(PROJECT_ROOT, "web", "dist", "index.html")

# A bun already on PATH at the pinned version, else the pinned release in bin/.
function ensure_bun()::String
    pin  = PINS["toolchain"]["bun"]
    want = string(pin["version"])
    for cand in (joinpath(BIN_DIR, "bun"), something(Sys.which("bun"), ""))
        (isempty(cand) || !isfile(cand)) && continue
        have = try strip(read(`$cand --version`, String)) catch; "" end
        have == want && return cand
    end
    rec = get(pin["archives"], PLATFORM, nothing)
    rec === nothing && error("No bun archive is pinned for $PLATFORM in $VERSIONS_FILE.")
    @info "Installing bun $want (pinned)..."
    zip = joinpath(BIN_DIR, "bun_download.zip")
    try
        download_verified(rec["url"], zip, get(rec, "sha256", nothing))
        return extract_binary(zip, "bun")
    finally
        rm(zip; force=true)
    end
end

# Stale when any frontend source, config or the lockfile is newer than the build.
function frontend_stale()::Bool
    isfile(DIST_INDEX) || return true
    built = mtime(DIST_INDEX)
    inputs = [joinpath(FRONTEND_DIR, f) for f in ("index.html", "package.json", "bun.lock", "vite.config.ts", "tsconfig.json")]
    for dir in ("src", "public"), (root, _, files) in walkdir(joinpath(FRONTEND_DIR, dir))
        append!(inputs, joinpath.(root, files))
    end
    any(p -> isfile(p) && mtime(p) > built, inputs)
end

function build_frontend(bun::String)
    cd(FRONTEND_DIR) do
        run(`$bun install --frozen-lockfile --ignore-scripts`)
        run(`$bun run build`)
    end
end

## Sysimage creation
const SYSIMAGE_EXT  = Sys.isapple() ? ".dylib" : ".so"
const SYSIMAGE_PATH = joinpath(PROJECT_ROOT, "MetaManifold$(SYSIMAGE_EXT)")
const PRECOMPILE_EXEC_PATH = joinpath(PROJECT_ROOT, "precompile_exec.jl")

function build_sysimage()
    # PackageCompiler is a project dependency, installed by Pkg.instantiate().
    @eval using PackageCompiler

    @info "Compiling sysimage - this may take several minutes...\n  Package: MetaManifold\n  Output:  $SYSIMAGE_PATH"

    kw = isfile(PRECOMPILE_EXEC_PATH) ?
        (; precompile_execution_file=[PRECOMPILE_EXEC_PATH]) : (;)

    PackageCompiler.create_sysimage(
        [:MetaManifold];
        sysimage_path = SYSIMAGE_PATH,
        project       = PROJECT_ROOT,
        kw...
    )

    @info "Sysimage written to $SYSIMAGE_PATH"
    println()
    println("Start the server with:")
    println("  bash start.sh")
end

## Preflight checks
"""Return true if the current process is running as root."""
has_root() = try; ccall(:geteuid, Cuint, ()) == 0; catch; false; end

"""Return true if sudo is available and can run non-interactively."""
function has_passwordless_sudo()::Bool
    Sys.which("sudo") === nothing && return false
    try
        run(pipeline(`sudo -n true`; stdout=devnull, stderr=devnull))
        return true
    catch
        return false
    end
end

function package_install_cmd(
    pkgs::Vector{String};
    brew_pkg::Union{String,Nothing}=nothing,
    pacman_pkg::Union{String,Nothing}=nothing,
    zypper_pkg::Union{String,Nothing}=nothing,
    linux_manager::Union{Symbol,Nothing}=nothing,
)
    if OS_TYPE == "macos"
        if Sys.which("brew") !== nothing
            mac_pkgs = brew_pkg === nothing ? pkgs : [brew_pkg]
            return Cmd(vcat(["brew", "install"], mac_pkgs))
        end
        return nothing
    end

    prefix = if has_root()
        String[]
    elseif has_passwordless_sudo()
        ["sudo"]
    else
        return nothing
    end
    manager = linux_manager
    if manager === nothing
        manager = Sys.which("apt") !== nothing ? :apt :
                  Sys.which("dnf") !== nothing ? :dnf :
                  Sys.which("pacman") !== nothing ? :pacman :
                  Sys.which("zypper") !== nothing ? :zypper :
                  nothing
    end
    manager === nothing && return nothing

    if manager == :apt
        return Cmd(vcat(prefix, ["apt", "install", "-y"], pkgs))
    elseif manager == :dnf
        return Cmd(vcat(prefix, ["dnf", "install", "-y"], pkgs))
    elseif manager == :pacman
        pacman_pkgs = pacman_pkg === nothing ? pkgs : [pacman_pkg]
        return Cmd(vcat(prefix, ["pacman", "-S", "--needed", "--noconfirm"], pacman_pkgs))
    elseif manager == :zypper
        zypper_pkgs = zypper_pkg === nothing ? pkgs : [zypper_pkg]
        return Cmd(vcat(prefix, ["zypper", "install", "-y"], zypper_pkgs))
    end

    nothing
end

## Main
function main()
    if UPDATE_MODE || MODIFY_MODE
        println()
        println("+-------------------------------------------+")
        println("|  MetaManifold - Install Script            |")
        println("+-------------------------------------------+")
        UPDATE_MODE && println("  Mode: UPDATE")
        MODIFY_MODE && println("  Mode: MODIFY")
        println()
    end

    # State the pins up front. An installer that says nothing about the versions it
    # is about to fetch leaves the operator no way to notice a wrong one.
    println("  Pinned tool versions ($VERSIONS_FILE):")
    for key in sort(collect(keys(PINS["tools"])))
        println("    $(rpad(bin_name(key), 12)) $(pinned_version(key))")
    end
    println()

    config = load_tools_config()
    resolved = Dict{String,Any}()

    # Resolve one tool and record how it landed, so the final summary can report
    # anything left unresolved instead of it surfacing later as a broken stage.
    function resolve_and_record(key, label, install_fn)
        p = resolve_tool(key, label, config, install_fn)
        resolved[key] = Dict("path" => p)
        p === nothing ?
            record!(label, STEP_ACTION, "unresolved - set a path under \"$key\" in config/tools.yml") :
            record!(label, STEP_OK, String(p))
        p
    end

    # Ensure pipx/pip is available before resolving Python-based tools
    ensure_pipx()

    resolve_and_record("cutadapt",   "cutadapt",   () -> install_python_tool("cutadapt"))
    # FastQC needs a Java runtime, so both are checked.
    fastqc_path = resolve_and_record("fastqc", "FastQC", () -> download_fastqc())
    fastqc_path === nothing || check_fastqc_runtime(fastqc_path)
    resolve_and_record("multiqc",    "MultiQC",    () -> install_python_tool("multiqc"))
    resolve_and_record("vsearch",    "vsearch",    () -> download_vsearch())
    resolve_and_record("cd_hit_est", "cd-hit-est", () -> download_cdhit())
    resolve_and_record("swarm",      "swarm",      () -> download_swarm())

    # Phylogenetic placement. MAFFT, IQ-TREE and RAxML may run on the bioserver
    # instead (pipeline.yml remote.stages); trimAl and gappa always run here.
    resolve_and_record("mafft",      "MAFFT",      () -> download_mafft())
    resolve_and_record("trimal",     "trimAl",     () -> download_trimal())
    resolve_and_record("iqtree",     "IQ-TREE",    () -> download_iqtree())
    resolve_and_record("raxml",      "RAxML",      () -> download_raxml())
    resolve_and_record("gappa",      "gappa",      () -> download_gappa())

    # Frontend. Built on every update, and on install whenever web/dist is
    # missing or older than the sources.
    println()
    println("  --- Frontend -------------------------------------------------------")
    try
        bun = ensure_bun()
        if UPDATE_MODE || frontend_stale()
            build_frontend(bun)
            record!("Frontend", STEP_OK, "built into web/dist with $bun")
        else
            record!("Frontend", STEP_OK, "web/dist is up to date")
        end
    catch e
        @error "Frontend build failed: $e"
        record!("Frontend", STEP_FAILED, "rerun install.sh, or build by hand: cd frontend && bun install && bun run build")
    end

    # R packages. Reproduced from renv.lock with renv::restore(), always, on a
    # plain install too - the DADA2 and NMDS/PERMANOVA stages need them. Writes
    # only to the project-local renv/library; never needs root.
    println()
    println("  --- R packages ----------------------------------------------------")
    if Sys.which("Rscript") === nothing
        @warn "Rscript not on PATH - skipping R package setup."
        record!("R runtime",  STEP_SKIPPED, "install R >= 4.0, then re-run install.sh")
        record!("R packages", STEP_SKIPPED, "blocked on the R runtime above")
    else
        rver = try
            strip(read(`Rscript -e "cat(as.character(getRversion()))"`, String))
        catch
            "unknown"
        end
        rpin = string(get(get(get(PINS, "runtimes", Dict()), "r", Dict()), "version", ""))
        if !isempty(rpin) && rver != "unknown" && !startswith(rver, rpin)
            record!("R runtime", STEP_OK, "R $rver (renv.lock pins $rpin; restore may rebuild from source)")
        else
            record!("R runtime", STEP_OK, "R $rver")
        end

        missing_hdrs = probe_missing_r_headers()
        if !isempty(missing_hdrs) && (has_root() || has_passwordless_sudo())
            # We already hold the privilege - use it, but never prompt for it.
            install_r_sysdeps()
            missing_hdrs = probe_missing_r_headers()
        end

        if !isempty(missing_hdrs)
            # renv::restore() from source would compile for minutes and then die on
            # the first package that needs one of these. Don't start it: report the
            # one command that unblocks it and stop here.
            hint = r_sysdep_hint(missing_hdrs)
            println("  Missing build dependencies - not running renv::restore().")
            println("    $hint")
            record!("R system headers", STEP_ACTION, hint)
            record!("R packages", STEP_ACTION,
                "renv::restore() skipped until the headers above are installed; then re-run install.sh")
        else
            rebuild = UPDATE_MODE && prompt_yn("  Force-rebuild every R package from source?", false)
            rstatus, rbad = setup_r_packages(; rebuild)
            rstatus == :ok ?
                record!("R packages", STEP_OK, "renv/library reproduced from renv.lock") :
                record!("R packages", STEP_ACTION,
                    "still unusable: $rbad - see the build errors above, then re-run: " *
                    "Rscript -e 'renv::restore(prompt = FALSE)'")
        end
    end

    # Write config
    write_tools_config(resolved)

    # Create a data directory
    mkpath("data")

    # Sysimage
    if SYSIMAGE_MODE
        println()
        println("  --- Julia sysimage -------------------------------------------------")
        try
            build_sysimage()
            record!("Julia sysimage", STEP_OK, SYSIMAGE_PATH)
        catch e
            @error "Sysimage build failed: $e"
            record!("Julia sysimage", STEP_FAILED, "rerun: julia --project=. install.jl --sysimage")
        end
    else
        record!("Julia sysimage", STEP_OK, "not requested (pass --sysimage to build one)")
    end

    println()
    println("Installation complete.")
    println()
    println("To adjust any tool locations, edit:")
    println("  $TOOLS_CONFIG")
    println("Stages that run on a server are set in the remote block of pipeline.yml.")

    print_summary()
end

RUNNING_AS_SCRIPT && main()
