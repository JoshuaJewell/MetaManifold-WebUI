#!/usr/bin/env julia
#
# Dependency installer for MetabarcodingPipeline
#
# Installs Julia deps, checks/downloads external CLI tools, and installs
# required R packages. Writes resolved tool paths to config/tools.yml.
#
# Usage via install.sh, or directly by:
#   julia --project=. install.jl [--update]
#
# Options:
#   --update    Re-check tool versions and update managed binaries in bin/

using Pkg
@info "Installing Julia package dependencies..."
Pkg.instantiate()

using YAML
import Downloads

## Instantiate
const UPDATE_MODE   = "--update" in ARGS
const PROJECT_ROOT  = @__DIR__
const BIN_DIR       = joinpath(PROJECT_ROOT, "bin")
const CONFIG_DIR    = joinpath(PROJECT_ROOT, "config")
const TOOLS_CONFIG  = joinpath(CONFIG_DIR, "tools.yml")

const OS_TYPE = Sys.islinux() ? "linux" :
                Sys.isapple() ? "macos" :
                error("Unsupported OS. Only Linux and macOS are supported.")

const ARCH_STR = Sys.ARCH == :x86_64  ? "x86_64"  :
                 Sys.ARCH == :aarch64 ? "aarch64" :
                 string(Sys.ARCH)

# Canonical binary name for tools whose config key differs from the binary name.
const BINARY_NAMES = Dict("cd_hit_est" => "cd-hit-est")
bin_name(key::String) = get(BINARY_NAMES, key, key)

mkpath(BIN_DIR)

## Config loading / saving
function load_tools_config()::Dict{String,Any}
    isfile(TOOLS_CONFIG) || return Dict{String,Any}()
    data = YAML.load_file(TOOLS_CONFIG)
    data isa Dict ? data : Dict{String,Any}()
end

function write_tools_config(config::Dict)
    mkpath(CONFIG_DIR)
    open(TOOLS_CONFIG, "w") do io
        for key in sort(collect(keys(config)))
            val = config[key]
            path = val isa Dict ? get(val, "path", nothing) : val
            println(io, "$key:")
            if path === nothing
                println(io, "  path: null")
            else
                # Quote the path to handle spaces and special characters
                println(io, "  path: \"$path\"")
            end
        end
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
# Fetch URL content as a String. Returns nothing on failure.
function fetch_string(url::String)::Union{String,Nothing}
    try
        buf = IOBuffer()
        Downloads.download(url, buf; headers=["User-Agent" => "MetabarcodingPipeline-installer"])
        String(take!(buf))
    catch e
        @warn "Could not fetch $url: $e"
        nothing
    end
end

# Return the first capture group of `pattern` in `s`, or nothing.
function first_match(pattern::Regex, s::String)::Union{String,Nothing}
    m = match(pattern, s)
    m === nothing ? nothing : m[1]
end

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

# Tool download functions

function download_vsearch()::String
    @info "Fetching latest vsearch release info from GitHub..."
    json = fetch_string("https://api.github.com/repos/torognes/vsearch/releases/latest")
    json === nothing && error("Cannot reach GitHub API. Check your internet connection.")

    # Find a download URL matching OS and architecture
    pattern = Regex(
        "\"browser_download_url\":\\s*\"(https://[^\"]*vsearch[^\"]*$(OS_TYPE)[^\"]*$(ARCH_STR)[^\"]*.tar.gz)\""
    )
    url = first_match(pattern, json)

    # Fallback to any tar.gz for this OS
    if url === nothing
        url = first_match(
            Regex("\"browser_download_url\":\\s*\"(https://[^\"]*vsearch[^\"]*$(OS_TYPE)[^\"]*.tar.gz)\""),
            json
        )
    end

    url === nothing && error(
        "No vsearch binary found for $OS_TYPE/$ARCH_STR in the latest release.\n" *
        "Check https://github.com/torognes/vsearch/releases for available builds."
    )

    tarball = joinpath(BIN_DIR, "vsearch_download.tar.gz")
    download_to(url, tarball)
    run(`tar -xzf $tarball -C $BIN_DIR`)
    rm(tarball)

    bin = find_file_in_dir(BIN_DIR, "vsearch")
    bin === nothing && error("vsearch binary not found after extraction.")

    dest = joinpath(BIN_DIR, "vsearch")
    bin != dest && mv(bin, dest; force=true)
    chmod(dest, 0o755)
    dest
end

function download_fastqc()::String
    # Pin to a known-good version; update periodically.
    version = "0.12.1"
    url = "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v$(version).zip"

    zipfile = joinpath(BIN_DIR, "fastqc.zip")
    download_to(url, zipfile)
    run(`unzip -q -o $zipfile -d $BIN_DIR`)
    rm(zipfile)

    bin = joinpath(BIN_DIR, "FastQC", "fastqc")
    isfile(bin) || error("fastqc not found after extraction. Expected at $bin")
    chmod(bin, 0o755)
    bin
end

function download_cdhit()::String
    # System package manager (no compilation needed)
    if OS_TYPE == "linux" && Sys.which("apt-get") !== nothing
        @info "Trying apt-get install cd-hit..."
        try
            run(`sudo apt-get install -y cd-hit`)
            found = Sys.which("cd-hit-est")
            found !== nothing && return found
        catch
            @warn "apt-get install failed (no sudo?), falling back to binary download."
        end
    elseif OS_TYPE == "macos" && Sys.which("brew") !== nothing
        @info "Trying brew install cd-hit..."
        try
            run(`brew install cd-hit`)
            found = Sys.which("cd-hit-est")
            found !== nothing && return found
        catch
            @warn "brew install failed, falling back to binary download."
        end
    end

    # Fallback to download precompiled binary.
    if OS_TYPE == "linux"
        url = "https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228-Linux.tar.gz"
    else
        error(
            "No managed cd-hit download available for macOS.\n" *
            "Install via:  brew install cd-hit\n" *
            "Then enter the path to cd-hit-est when prompted."
        )
    end

    @info "Downloading cd-hit v4.8.1 precompiled binary..."
    tarball = joinpath(BIN_DIR, "cdhit_download.tar.gz")
    download_to(url, tarball)
    run(`tar -xzf $tarball -C $BIN_DIR`)
    rm(tarball)

    bin = find_file_in_dir(BIN_DIR, "cd-hit-est")
    bin === nothing && error("cd-hit-est binary not found after extraction.")

    dest = joinpath(BIN_DIR, "cd-hit-est")
    bin != dest && mv(bin, dest; force=true)
    chmod(dest, 0o755)
    dest
end

function install_python_tool(name::String)::String
    # pipx (recommended on PEP 668 / Debian-managed systems)
    if Sys.which("pipx") !== nothing
        @info "Installing $name via pipx..."
        run(`pipx install $name`)
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

    @info "Installing $name via pip..."
    success = try
        run(`$pip_cmd install --user $name`)
        true
    catch
        false
    end

    if !success
        # Fallback 2 to PEP 668: externally-managed environment, try --break-system-packages
        @warn "pip --user blocked by system policy. Retrying with --break-system-packages..."
        run(`$pip_cmd install --user --break-system-packages $name`)
    end

    # Locate the installed binary
    found = find_in_path(name)
    found !== nothing && return found

    local_bin = joinpath(homedir(), ".local", "bin", name)
    isfile(local_bin) && return local_bin

    @warn "Installed $name but could not locate the binary. Ensure ~/.local/bin is in PATH."
    name
end

function install_r_packages(packages::Vector{String})
    pkgs_r = join(["\"$p\"" for p in packages], ", ")
    snippet = """
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager", repos = "https://cloud.r-project.org")
        pkgs <- c($pkgs_r)
        missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
        if (length(missing) > 0) {
            message("Installing: ", paste(missing, collapse = ", "))
            BiocManager::install(missing, ask = FALSE)
        } else {
            message("All R packages already installed.")
        }
    """
    try
        run(`Rscript -e $snippet`)
        @info "R packages installed successfully."
    catch e
        @error "R package installation failed: $e"
    end
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
    println()
    println("  ─── $display_name ─────────────────────────────────────────")

    bin = bin_name(key)

    # Check existing config (skip in update mode for managed installs)
    existing = get(existing_config, key, nothing)
    existing_path = existing isa Dict ? get(existing, "path", nothing) : nothing

    if existing_path !== nothing && !UPDATE_MODE
        println("  Configured path: $existing_path")
        if check_tool(string(existing_path))
            prompt_yn("  Use this?") && return string(existing_path)
        else
            println("  Warning: configured path does not appear to be callable.")
        end
    end

    # Check PATH
    path_result = find_in_path(bin)
    if path_result !== nothing
        println("  Found in PATH:   $path_result")
        prompt_yn("  Use this?") && return path_result
    end

    # Build option list
    options = String[]
    if install_fn !== nothing
        push!(options, "Install/download automatically to bin/")
    end
    push!(options, "Enter a path manually  (local: /path/to/$bin  or  remote: user@host:/path/to/$bin)")
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
                println("  Falling back to manual entry.")
                choice = 2
            end
        end
        # After potential fallback, re-map remaining choices
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

## Main
function main()
    println()
    println("╔═══════════════════════════════════════════╗")
    println("║  MetabarcodingPipeline — Install Script   ║")
    println("╚═══════════════════════════════════════════╝")
    UPDATE_MODE && println("  Mode: UPDATE")
    println()

    config = load_tools_config()
    resolved = Dict{String,Any}()

    # cutadapt
    path = resolve_tool("cutadapt", "cutadapt", config,
        () -> install_python_tool("cutadapt"))
    resolved["cutadapt"] = Dict("path" => path)

    # fastqc
    path = resolve_tool("fastqc", "FastQC", config,
        () -> download_fastqc())
    resolved["fastqc"] = Dict("path" => path)

    # multiqc
    path = resolve_tool("multiqc", "MultiQC", config,
        () -> install_python_tool("multiqc"))
    resolved["multiqc"] = Dict("path" => path)

    # vsearch
    path = resolve_tool("vsearch", "vsearch", config,
        () -> download_vsearch())
    resolved["vsearch"] = Dict("path" => path)

    # cd-hit-est
    path = resolve_tool("cd_hit_est", "cd-hit-est", config,
        () -> download_cdhit())
    resolved["cd_hit_est"] = Dict("path" => path)

    # R packages
    println()
    println("  ─── R packages ─────────────────────────────────────────────")
    r_packages = ["dada2", "tidyverse"]
    if prompt_yn("  Install/check R packages (dada2, tidyverse)?")
        install_r_packages(r_packages)
    end

    # Write config
    write_tools_config(resolved)

    # Create a data directory
    mkpath("data")

    println()
    println("Installation complete.")
    println()
    println("To set remote SSH paths or adjust any tool locations, edit:")
    println("  $TOOLS_CONFIG")
    println()
    println("Example remote path (SSH):")
    println("  vsearch:")
    println("    path: \"user@bioserver:/home/user/software/vsearch\"")
end

main()