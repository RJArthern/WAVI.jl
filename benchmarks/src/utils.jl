using Dates
using JSON3
using MPI
using Printf
using Profile
using Logging
using LoggingExtras

const BENCHMARK_OUTPUT_DIR = normpath(@__DIR__, "..", "output")

"""
    BenchmarkResults

Summary of one timed driver run: wall time (local and MPI-max), allocations, GC,
peak RSS (local, MPI-max, and MPI-sum), and the sample interval used for the
resource time series.
"""
struct BenchmarkResults
    execution_time::Float64
    execution_time_max::Float64
    allocated_bytes::Int64
    gc_time::Float64
    allocations::Int64
    peak_rss_bytes::Int
    peak_rss_max_bytes::Int
    peak_rss_sum_bytes::Int
    sample_interval_s::Float64
    system_info::Dict{String, Any}
    timestamp::DateTime
end

"""
    write_resource_timeseries(path, samples)

Write resource samples to CSV.

Non-MPI (and single-rank MPI) columns:
`elapsed_s,cpu_cores_used,rss_max_bytes,rss_sum_bytes,cpu_rank0,rss_rank0`
(`rss_max` / `rss_sum` match the local process).

When running with MPI, all processors must call this function together. Rank 0
collects everyone's CPU and memory data and saves it into a single CSV file.
This file shows the total resource usage across the entire cluster, alongside
individual columns for each processor (e.g., `cpu_rank0`, `cpu_rank1`).
"""
function write_resource_timeseries(path::AbstractString, samples)
    if MPI.Initialized() && MPI.Comm_size(MPI.COMM_WORLD) > 1
        _write_resource_timeseries_mpi(path, samples)
    else
        _write_resource_timeseries_local(path, samples)
    end
    return nothing
end

function _write_resource_timeseries_local(path::AbstractString, samples)
    open(path, "w") do io
        println(io, "elapsed_s,cpu_cores_used,rss_max_bytes,rss_sum_bytes,cpu_rank0,rss_rank0")
        for s in samples
            println(
                io,
                "$(s.elapsed_s),$(s.cpu_cores_used),$(s.rss_bytes),$(s.rss_bytes),$(s.cpu_cores_used),$(s.rss_bytes)",
            )
        end
    end
end

function _write_resource_timeseries_mpi(path::AbstractString, samples)
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)
    n = Int(MPI.Allreduce(length(samples), MPI.MIN, comm))

    header_parts = String["elapsed_s", "cpu_cores_used", "rss_max_bytes", "rss_sum_bytes"]
    append!(header_parts, ["cpu_rank$r" for r in 0:(nranks - 1)])
    append!(header_parts, ["rss_rank$r" for r in 0:(nranks - 1)])
    header = join(header_parts, ',')

    if n == 0
        if rank == 0
            open(path, "w") do io
                println(io, header)
            end
        end
        return nothing
    end

    elapsed = Vector{Float64}(undef, n)
    cpu = Vector{Float64}(undef, n)
    rss = Vector{Int64}(undef, n)
    @inbounds for i in 1:n
        elapsed[i] = samples[i].elapsed_s
        cpu[i] = samples[i].cpu_cores_used
        rss[i] = Int64(samples[i].rss_bytes)
    end

    cpu_flat = MPI.Gather(cpu, 0, comm)
    rss_flat = MPI.Gather(rss, 0, comm)
    rank == 0 || return nothing

    cpu_mat = reshape(cpu_flat, n, nranks)
    rss_mat = reshape(rss_flat, n, nranks)

    open(path, "w") do io
        println(io, header)
        @inbounds for i in 1:n
            cpu_row = @view cpu_mat[i, :]
            rss_row = @view rss_mat[i, :]
            cpu_sum = sum(cpu_row)
            rss_max = maximum(rss_row)
            rss_sum = sum(rss_row)
            print(io, elapsed[i], ',', cpu_sum, ',', rss_max, ',', rss_sum)
            for r in 1:nranks
                print(io, ',', cpu_mat[i, r])
            end
            for r in 1:nranks
                print(io, ',', rss_mat[i, r])
            end
            println(io)
        end
    end
    return nothing
end

"""
    benchmark_main(id, model, model_args, variables_to_plot, rank=0; kwargs...)

Create an output folder, optionally warm up, time `model`, save JSON/CSV, and
optionally plot NetCDF fields. Only rank 0 writes files and plots.
"""
function benchmark_main(id::String,
                        model::Function,
                        model_args::Dict,
                        variables_to_plot::Vector{String},
                        rank::Int = 0;
                        metadata::Dict{String, Any} = Dict{String, Any}(),
                        sample_interval::Float64 = 0.25,
                        no_plots::Bool = false,
                        warmup::Bool = false)

    driver_name = get(metadata, "driver", "unknown_driver")
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    if MPI.Initialized()
        # Ensure all ranks use the exact same timestamp string to prevent
        # split output directories if rank clocks are slightly desynced
        timestamp_chars = collect(timestamp)
        MPI.Bcast!(timestamp_chars, 0, MPI.COMM_WORLD)
        timestamp = String(timestamp_chars)
    end

    output_dir = joinpath(
        BENCHMARK_OUTPUT_DIR,
        driver_name,
        "benchmark_$(id)_$(timestamp)",
    )
    model_args[:folder] = output_dir
    mkpath(output_dir)

    old_logger = global_logger()
    log_io = nothing
    if rank == 0
        log_io = open(joinpath(output_dir, "run.log"), "w")
        global_logger(TeeLogger(old_logger, SimpleLogger(log_io)))
        @async begin
            while isopen(log_io)
                flush(log_io)
                sleep(1)
            end
        end
    end

    try

    # Copy the driver adaptor file to the output directory for reproducibility
    driver_path = joinpath(normpath(@__DIR__, "..", "drivers"), "$(driver_name).jl")
    if isfile(driver_path) && rank == 0
        cp(driver_path, joinpath(output_dir, "$(driver_name).jl"); force=true)
    end

    if rank == 0
        @info "JULIA WAVI BENCHMARK SUITE - $(id)"
        @info "Start time: $(now())"
        if haskey(metadata, "driver")
            @info "Driver: $(metadata["driver"])"
        end
        if haskey(metadata, "command") && !isnothing(metadata["command"])
            @info "Command: $(metadata["command"])"
        end
        @info "Executing model with $(model_args)"
    end

    if warmup
        rank == 0 && @info "Warm-up run (not timed)..."
        Base.invokelatest(model; model_args...)
    end

    # All ranks pass the same path so MPI sample gather stays collective; only
    # rank 0 writes the CSV inside `write_resource_timeseries`.
    timeseries_path = joinpath(output_dir, "resource_timeseries.csv")
    result, benchmark_results = monitor_resources(
        model;
        sample_interval = sample_interval,
        timeseries_path = timeseries_path,
        model_args...,
    )

    if rank == 0
        # Display results
        @info "Execution time: $(@sprintf("%.3f", benchmark_results.execution_time)) seconds"
        if benchmark_results.execution_time_max > benchmark_results.execution_time
            @info "Execution time (MPI max): $(@sprintf("%.3f", benchmark_results.execution_time_max)) seconds"
        end
        @info "Allocated: $(@sprintf("%.2f", benchmark_results.allocated_bytes / 1024^2)) MB"
        @info "Peak RSS: $(@sprintf("%.2f", benchmark_results.peak_rss_bytes / 1024^2)) MB"
        if benchmark_results.peak_rss_max_bytes > benchmark_results.peak_rss_bytes ||
           benchmark_results.peak_rss_sum_bytes > benchmark_results.peak_rss_bytes
            @info "Peak RSS (MPI max): $(@sprintf("%.2f", benchmark_results.peak_rss_max_bytes / 1024^2)) MB"
            @info "Peak RSS (MPI sum): $(@sprintf("%.2f", benchmark_results.peak_rss_sum_bytes / 1024^2)) MB"
        end
        @info "GC time: $(@sprintf("%.3f", benchmark_results.gc_time)) seconds"
        @info "Allocations: $(benchmark_results.allocations)"

        if hasproperty(result, :setup_time)
            @info "Setup time: $(@sprintf("%.3f", result.setup_time)) seconds"
            metadata["setup_time_seconds"] = result.setup_time
        end
        if hasproperty(result, :solve_time)
            @info "Solve time: $(@sprintf("%.3f", result.solve_time)) seconds"
            metadata["solve_time_seconds"] = result.solve_time
        end
        if hasproperty(result, :setup_time) && hasproperty(result, :solve_time)
            comp_time = max(0.0, benchmark_results.execution_time - (result.setup_time + result.solve_time))
            @info "Compilation & Overhead time: $(@sprintf("%.3f", comp_time)) seconds"
            metadata["compilation_and_overhead_time_seconds"] = comp_time
        end

        benchmark_file = joinpath(output_dir, "benchmark_results.json")
        save_benchmark_results(benchmark_results, benchmark_file; metadata = metadata)

        netcdf_output = joinpath(output_dir, "outfile.nc")
        if !no_plots && isfile(netcdf_output)
            @info "Creating visualisations..."
            plot_multiple_variables(netcdf_output, variables_to_plot, joinpath(output_dir, "plots"))
        elseif no_plots
            @info "Skipping plots (--no-plots)."
        else
            @info "Warning: NetCDF output file '$netcdf_output' not found. Skipping visualisation."
        end

        @info "End time: $(now())"
        @info "All output saved to: $output_dir"
    end

    return result, benchmark_results
    finally
        if rank == 0 && log_io !== nothing
            global_logger(old_logger)
            close(log_io)
        end
    end
end

"""
    monitor_resources(func, args...; sample_interval=0.25, timeseries_path=nothing, kwargs...)

Time `func(args...; kwargs...)` while sampling RSS and CPU.

The resource monitor starts before `@timed` and stops afterwards, so helper
allocations are not counted in `allocated_bytes`. If `timeseries_path` is set,
samples are written there as CSV.

When MPI is initialised, wall time is reduced with `MPI.MAX` so every rank
shares the slowest-rank time as `execution_time_max`. Peak RSS is reduced with
`MPI.MAX` and `MPI.SUM` for cluster-wide memory figures. The resource CSV is
likewise gathered across ranks (see `write_resource_timeseries`).
"""
function monitor_resources(func, args...;
                           sample_interval::Float64 = 0.25,
                           timeseries_path::Union{Nothing, AbstractString} = nothing,
                           kwargs...)
    Profile.clear()

    mon = start_monitor!(sample_interval)
    local result
    try
        @profile begin
            # invokelatest: driver functions are included at runtime (Fixes world age error)
            result = @timed Base.invokelatest(func, args...; kwargs...)
        end
    finally
        stop_monitor!(mon)
    end

    if timeseries_path !== nothing
        write_resource_timeseries(timeseries_path, mon.samples)
        if !MPI.Initialized() || MPI.Comm_rank(MPI.COMM_WORLD) == 0
            @info "Resource time series saved to: $timeseries_path"
        end
    end

    t_local = result.time
    rss_local = mon.peak_rss_bytes
    if MPI.Initialized()
        t_max = MPI.Allreduce(t_local, MPI.MAX, MPI.COMM_WORLD)
        rss_max = Int(MPI.Allreduce(rss_local, MPI.MAX, MPI.COMM_WORLD))
        rss_sum = Int(MPI.Allreduce(rss_local, MPI.SUM, MPI.COMM_WORLD))
    else
        t_max = t_local
        rss_max = rss_local
        rss_sum = rss_local
    end

    benchmark_result = BenchmarkResults(
        t_local,
        t_max,
        Int64(result.bytes),
        result.gctime,
        Int64(result.gcstats.allocd),
        rss_local,
        rss_max,
        rss_sum,
        sample_interval,
        get_system_info(),
        now(),
    )

    return result.value, benchmark_result
end

function get_system_info()
    return Dict(
        "julia_version" => string(VERSION),
        "cpu_threads" => Sys.CPU_THREADS,
        "total_memory" => Sys.total_memory(),
        "hostname" => gethostname(),
        "os" => string(Sys.KERNEL),
    )
end

"""
    slurm_metadata() -> Dict{String, Any}

Reads SLURM env vars if they are set.

When run inside a SLURM allocation, this copies scheduler env vars into the
results JSON with a `slurm_` prefix, e.g. `slurm_job_id`, `slurm_ntasks`.
Outside SLURM, this returns an empty dict so local runs stay unchanged.
"""
function slurm_metadata()
    mapping = (
        "SLURM_JOB_ID" => "slurm_job_id",
        "SLURM_NNODES" => "slurm_nnodes",
        "SLURM_NTASKS" => "slurm_ntasks",
        "SLURM_NTASKS_PER_NODE" => "slurm_ntasks_per_node",
        "SLURM_CPUS_PER_TASK" => "slurm_cpus_per_task",
        "SLURM_JOB_PARTITION" => "slurm_job_partition",
        "SLURM_JOB_QOS" => "slurm_job_qos",
        "SLURM_JOB_ACCOUNT" => "slurm_job_account",
    )

    out = Dict{String, Any}()
    for (env_name, key) in mapping
        haskey(ENV, env_name) || continue
        raw = ENV[env_name]
        parsed = tryparse(Int, raw)
        out[key] = parsed === nothing ? raw : parsed
    end
    return out
end

"""
    save_benchmark_results(results, filename; metadata=Dict())

Write timing and memory summary fields to JSON, including `metadata`.
"""
function save_benchmark_results(results::BenchmarkResults, filename::String;
                                metadata::Dict{String, Any} = Dict{String, Any}())
    output_data = Dict(
        "timestamp" => string(results.timestamp),
        "execution_time_seconds" => results.execution_time,
        "execution_time_max_seconds" => results.execution_time_max,
        "allocated_bytes" => results.allocated_bytes,
        "peak_rss_bytes" => results.peak_rss_bytes,
        "peak_rss_max_rank_bytes" => results.peak_rss_max_bytes,
        "peak_rss_sum_ranks_bytes" => results.peak_rss_sum_bytes,
        "sample_interval_s" => results.sample_interval_s,
        "gc_time_seconds" => results.gc_time,
        "allocations" => results.allocations,
        "system_info" => results.system_info,
        "metadata" => metadata,
    )

    open(filename, "w") do io
        JSON3.pretty(io, output_data)
    end

    @info "Benchmark results saved to: $filename"
end
