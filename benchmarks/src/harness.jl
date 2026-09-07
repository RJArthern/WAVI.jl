# # WAVI.jl benchmark harness

using LinearAlgebra
using MPI
using WAVI

const _VALID_MODES = (:basic, :threaded, :mpi)

# # Configuration

"""
    BenchmarkOptions

Configuration for a benchmark run.

# Fields

- `mode`: execution mode: `:basic`, `:threaded`, or `:mpi`
- `driver`: registered adaptor name (e.g. `"mismip_plus"`)
- `ngridsx`, `ngridsy`, `overlap`, `niterations`: ThreadedSpec / Schwarz parameters
- `px`, `py`: MPI process grid dimensions (`px == 0` means use `Comm_size`; `py` defaults to `1`, i.e. an `N×1` layout)
- `sample_interval`: seconds between resource samples
- `no_plots`: skip NetCDF plots if true
- `warmup`: untimed run before the measured benchmark
"""
Base.@kwdef struct BenchmarkOptions
    mode::Symbol
    driver::String
    ngridsx::Int = 2
    ngridsy::Int = 2
    overlap::Int = 2
    niterations::Int = 2
    px::Int = 0
    py::Int = 1
    sample_interval::Float64 = 0.25
    no_plots::Bool = false
    warmup::Bool = false
    tag::String = ""
    output_group::String = ""
end

# Convenience constructor accepting a string mode (from the Comonicon CLI),
# coerced to the lower-case Symbol the dispatch expects.
BenchmarkOptions(mode::AbstractString; kwargs...) =
    BenchmarkOptions(; mode = Symbol(lowercase(mode)), kwargs...)

# Set from benchmarks/run.jl before the CLI runs; logged at the start of run_benchmark.
const BENCHMARK_COMMAND = Ref{Union{String, Nothing}}(nothing)

"""
    set_benchmark_command!(cmd::AbstractString)

Store the exact command-line invocation used to launch the benchmark.
This is logged to the console at the start of `run_benchmark` and is intended
to be saved alongside benchmark metadata for reproducibility.
"""
set_benchmark_command!(cmd::AbstractString) = (BENCHMARK_COMMAND[] = String(cmd))

# # BLAS / threading

"""
    pin_blas_threads!()

Pin BLAS to exactly one thread per process by setting OpenBLAS/MKL environment variables.
This is critical to prevent severe thread oversubscription during MPI or ThreadedSpec runs.
"""
function pin_blas_threads!()
    !haskey(ENV, "OPENBLAS_NUM_THREADS") && (ENV["OPENBLAS_NUM_THREADS"] = "1")
    !haskey(ENV, "MKL_NUM_THREADS") && (ENV["MKL_NUM_THREADS"] = "1")
    BLAS.set_num_threads(1)
end

"""
    benchmark_metadata(opts::BenchmarkOptions; mpi_world_size=nothing)

Extract benchmark details (mode, driver, julia threads, and BLAS threads)
into a dict for serialisation with final benchmark JSON results output.
"""
function benchmark_metadata(opts::BenchmarkOptions; mpi_world_size::Union{Nothing, Int} = nothing)
    commit_hash = ""
    try
        commit_hash = strip(read(`git rev-parse HEAD`, String))
    catch
        commit_hash = "unknown"
    end

    return merge!(Dict{String, Any}(
        "mode" => string(opts.mode),
        "driver" => opts.driver,
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "sample_interval_s" => opts.sample_interval,
        "reference_cores" => reference_cores(opts; mpi_world_size = mpi_world_size),
        "command" => BENCHMARK_COMMAND[],
        "tag" => opts.tag,
        "output_group" => opts.output_group,
        "git_commit" => commit_hash,
    ), slurm_metadata())
end

"""
    reference_cores(opts; mpi_world_size=nothing) -> Int

How many cores this benchmark is meant to use, for normalising CPU samples
(`cpu_fraction = cpu_cores_used / reference_cores`).

- `:basic`: `Threads.nthreads()`
- `:threaded`: `ngridsx * ngridsy`
- `:mpi`: `mpi_world_size` if given, else `SLURM_NTASKS`, else `1`
"""
function reference_cores(opts::BenchmarkOptions; mpi_world_size::Union{Nothing, Int} = nothing)
    if opts.mode == :threaded
        return opts.ngridsx * opts.ngridsy
    elseif opts.mode == :mpi
        return something(mpi_world_size, 1)
    else
        return Threads.nthreads()
    end
end

# # Execution dispatch

"""
    run_benchmark(opts::BenchmarkOptions)

Load the requested driver adaptor, build appropriate WAVI spec, and run
`benchmark_main` for monitoring.
"""
function run_benchmark(opts::BenchmarkOptions)
    opts.mode in _VALID_MODES || error("Unknown mode '$(opts.mode)'. Use basic, threaded, or mpi.")

    pin_blas_threads!()

    # Initialise MPI up front so the entire setup is covered by finalisation.
    opts.mode == :mpi && !MPI.Initialized() && MPI.Init()

    rank = opts.mode == :mpi ? MPI.Comm_rank(MPI.COMM_WORLD) : 0
    mpi_world_size = opts.mode == :mpi ? MPI.Comm_size(MPI.COMM_WORLD) : nothing
    metadata = benchmark_metadata(opts; mpi_world_size = mpi_world_size)



    driver = load_driver(opts.driver)

    try
        run_id = ""
        spec_kwargs = Dict{Symbol, Any}()

        if opts.mode == :basic
            # Serial / Multi-threaded: BasicSpec setup
            run_id = "basic" * (Threads.nthreads() > 1 ? "_t$(Threads.nthreads())" : "")

        elseif opts.mode == :threaded
            # Shared memory: ThreadedSpec setup
            nx, ny, ov, ni = opts.ngridsx, opts.ngridsy, opts.overlap, opts.niterations
            spec_kwargs[:spec] = ThreadedSpec(; ngridsx = nx, ngridsy = ny, overlap = ov, niterations = ni)
            run_id = "threaded.$(nx)x$(ny)_o$(ov)_i$(ni)"

        elseif opts.mode == :mpi
            # Distributed: MPISpec setup
            comm = MPI.COMM_WORLD
            sz = MPI.Comm_size(comm)

            px = opts.px == 0 ? sz : opts.px
            py = opts.py
            px * py == sz || error("MPI process grid px×py ($(px)×$(py)) must equal world size ($(sz)).")


            grid = Base.invokelatest(driver.grid)
            # Narrow domains (e.g. MISMIP+ ny=10) need enough core cells after halo.
            # Prefer px×1 (default py=1); warn when py>1 leaves a thin y-patch.
            halo = 2
            if rank == 0 && py > 1
                core_y = div(grid.ny, py)
                if core_y <= 2 * halo
                    @warn "MPI process grid $(px)×$(py) leaves ~$(core_y) core cells in y " *
                          "(halo=$halo, ny=$(grid.ny)). Prefer a 1D layout (e.g. --px $sz --py 1) " *
                          "for PoU on narrow domains such as MISMIP+."
                end
            end
            spec = MPISpec(px, py, halo, grid; pou = true, niterations = opts.niterations)

            spec_kwargs[:grid] = grid
            spec_kwargs[:spec] = spec
            run_id = "mpi.$(px)x$(py)_sz$(sz)"
        end

        benchmark_main(run_id, driver.run, spec_kwargs, driver.plot_vars, rank;
                       metadata = metadata,
                       sample_interval = opts.sample_interval,
                       no_plots = opts.no_plots,
                       warmup = opts.warmup)
    finally
        opts.mode == :mpi && MPI.Initialized() && MPI.Finalize()
    end

    return nothing
end

"""
    log_wavelet_sizes(result)

Print coarse wavelet DOFs (`wu.n + wv.n`) against fine inner DOFs (`gu.ni + gv.ni`).
Used after a profile warm-up to decide whether assembling a coarse operator is cheap.
If `io` is given, the same line is written there as well.
"""
function log_wavelet_sizes(result; io::Union{IO, Nothing} = nothing)
    sim = if result isa NamedTuple && haskey(result, :simulation)
        result.simulation
    elseif hasproperty(result, :model)
        result
    else
        nothing
    end
    sim === nothing && return
    hasproperty(sim, :model) || return
    fields = sim.model.fields
    n_fine = fields.gu.ni + fields.gv.ni
    n_coarse = fields.wu.n[] + fields.wv.n[]
    msg = "Wavelet sizes: n_coarse=$(n_coarse) (wu.n=$(fields.wu.n[]), wv.n=$(fields.wv.n[])); " *
          "n_fine=$(n_fine) (gu.ni=$(fields.gu.ni), gv.ni=$(fields.gv.ni)); " *
          "ratio=$(round(n_coarse / max(n_fine, 1); digits = 4))"
    @info msg
    io !== nothing && println(io, msg)
    return nothing
end

"""
    run_profile(opts::BenchmarkOptions)

Warm up once, then `@profile` the driver. Rank 0 writes a flat Profile report
under `benchmarks/output/profile_<id>_<timestamp>/`.

For allocation profiling, launch Julia with `--track-allocation=user` instead
(slow; not part of this subcommand).
"""
function run_profile(opts::BenchmarkOptions)
    opts.mode in _VALID_MODES || error("Unknown mode '$(opts.mode)'. Use basic, threaded, or mpi.")

    pin_blas_threads!()
    opts.mode == :mpi && !MPI.Initialized() && MPI.Init()

    rank = opts.mode == :mpi ? MPI.Comm_rank(MPI.COMM_WORLD) : 0
    rank == 0 && !isnothing(BENCHMARK_COMMAND[]) && @info "Command: $(BENCHMARK_COMMAND[])"

    driver = load_driver(opts.driver)

    try
        spec_kwargs = Dict{Symbol, Any}()
        run_id = "basic"

        if opts.mode == :threaded
            nx, ny, ov, ni = opts.ngridsx, opts.ngridsy, opts.overlap, opts.niterations
            spec_kwargs[:spec] = ThreadedSpec(; ngridsx = nx, ngridsy = ny, overlap = ov, niterations = ni)
            run_id = "threaded.$(nx)x$(ny)_o$(ov)_i$(ni)"

        elseif opts.mode == :mpi
            sz = MPI.Comm_size(MPI.COMM_WORLD)
            px = opts.px == 0 ? sz : opts.px
            py = opts.py
            px * py == sz || error("MPI process grid px×py ($(px)×$(py)) must equal world size ($(sz)).")


            grid = Base.invokelatest(driver.grid)
            halo = 2
            if rank == 0 && py > 1
                core_y = div(grid.ny, py)
                if core_y <= 2 * halo
                    @warn "MPI process grid $(px)×$(py) leaves ~$(core_y) core cells in y " *
                          "(halo=$halo, ny=$(grid.ny)). Prefer a 1D layout (e.g. --px $sz --py 1) " *
                          "for PoU on narrow domains such as MISMIP+."
                end
            end
            spec_kwargs[:grid] = grid
            spec_kwargs[:spec] = MPISpec(px, py, halo, grid; pou = false, niterations = opts.niterations)
            run_id = "mpi.$(px)x$(py)_sz$(sz)"
        end

        output_dir = joinpath(
            BENCHMARK_OUTPUT_DIR,
            "profile_$(run_id)_$(Dates.format(Dates.now(), "yyyymmdd_HHMMSS"))",
        )
        mkpath(output_dir)
        spec_kwargs[:folder] = output_dir

        rank == 0 && @info "Profiling $(opts.driver) ($(opts.mode)): warm-up then @profile..."
        warmup_result = Base.invokelatest(driver.run; spec_kwargs...)
        if rank == 0
            open(joinpath(output_dir, "wavelet_sizes.txt"), "w") do io
                log_wavelet_sizes(warmup_result; io = io)
            end
        end
        # Default n=10^7 fills on the 381x381 ISMIP7 driver; 10^8 is about 800 MB and
        # covers a full 5-step run at -t 1 and a moderately threaded -t 16 run.
        Profile.init(n = 10^8, delay = 0.001)
        Profile.clear()
        @profile Base.invokelatest(driver.run; spec_kwargs...)

        if rank == 0
            profile_file = joinpath(output_dir, "profile_flat.txt")
            open(profile_file, "w") do io
                Profile.print(io; format = :flat)
            end
            @info "Flat profile saved to: $profile_file"
        end
    finally
        opts.mode == :mpi && MPI.Initialized() && MPI.Finalize()
    end

    return nothing
end
