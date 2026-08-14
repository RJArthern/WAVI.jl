# WAVIBenchmarks — benchmark harness for WAVI.jl (BasicSpec / ThreadedSpec / MPISpec).
#
# To install dependencies the first time, or on Project.toml changes, run:
#   cd benchmarks
#   julia --project=. -e 'using Pkg; Pkg.instantiate()'

module WAVIBenchmarks

using Comonicon

# Include shared harness components
include(joinpath(@__DIR__, "drivers.jl"))
include(joinpath(@__DIR__, "resource_monitor.jl"))
include(joinpath(@__DIR__, "utils.jl"))
include(joinpath(@__DIR__, "plotting.jl"))
include(joinpath(@__DIR__, "harness.jl"))

export BenchmarkOptions, run_benchmark, run_profile, set_benchmark_command!

"""
List registered benchmark driver adaptors discovered under `benchmarks/drivers/`.
This command checks the registry for available drivers that can be used 
with the benchmark harness and prints them to the console.
"""
@cast function drivers()
    names = available_drivers()
    println("Available drivers: ", isempty(names) ? "(none)" : join(names, ", "))
end

"""
Run a timed benchmark for the given execution mode and driver adaptor.

# Arguments

- `mode`: `basic`, `threaded`, or `mpi`
- `driver`: registered adaptor name (e.g. `mismip_plus`)

# Options

- `--niterations <n>`: [ThreadedSpec, MPISpec] Schwarz/PoU solver iterations (default: 2)

- `--ngridsx <n>`: [ThreadedSpec] x domain decomposition (default: 2)
- `--ngridsy <n>`: [ThreadedSpec] y domain decomposition (default: 2)
- `--overlap <n>`: [ThreadedSpec] Schwarz overlap cells (default: 2)

- `--px <n>`: [MPISpec] process grid x (`0` = MPI world size; default: 0 → `N×1` with `--py 1`)
- `--py <n>`: [MPISpec] process grid y (default: 1; prefer `1` on narrow domains like MISMIP+)

- `--sample-interval <s>`: resource sample period in seconds (default: 0.25)
- `--no-plots`: skip NetCDF plots
- `--warmup`: run once untimed before the measured run
"""
@cast function run(
    mode::String,
    driver::String;
    ngridsx::Int = 2,
    ngridsy::Int = 2,
    overlap::Int = 2,
    niterations::Int = 2,
    px::Int = 0,
    py::Int = 1,
    sample_interval::Float64 = 0.25,
    no_plots::Bool = false,
    warmup::Bool = false,
    tag::String = "",
)
    opts = BenchmarkOptions(
        mode;
        driver = driver,
        ngridsx = ngridsx,
        ngridsy = ngridsy,
        overlap = overlap,
        niterations = niterations,
        px = px,
        py = py,
        sample_interval = sample_interval,
        no_plots = no_plots,
        warmup = warmup,
        tag = tag,
    )
    run_benchmark(opts)
end

"""
On-demand CPU profile of a driver (separate from timed `run`).
Warm up once, then `@profile`, and write a flat report under `benchmarks/output/`.
For allocation profiling, launch Julia with `--track-allocation=user` instead
(slow; not handled by this subcommand).

# Arguments

- `driver`: registered adaptor name (e.g. `mismip_plus`)

# Options

- `--mode <mode>`: `basic` (default), `threaded`, or `mpi`
- `--niterations <n>`: [ThreadedSpec, MPISpec] Schwarz/PoU iterations (default: 2)
- `--ngridsx <n>`: [ThreadedSpec] x domain decomposition (default: 2)
- `--ngridsy <n>`: [ThreadedSpec] y domain decomposition (default: 2)
- `--overlap <n>`: [ThreadedSpec] Schwarz overlap cells (default: 2)
- `--px <n>`: [MPISpec] process grid x (`0` = MPI world size; default: 0)
- `--py <n>`: [MPISpec] process grid y (default: 1; prefer `1` on narrow domains)
"""
@cast function profile(
    driver::String;
    mode::String = "basic",
    ngridsx::Int = 2,
    ngridsy::Int = 2,
    overlap::Int = 2,
    niterations::Int = 2,
    px::Int = 0,
    py::Int = 1,
    tag::String = "",
)
    opts = BenchmarkOptions(
        mode;
        driver = driver,
        ngridsx = ngridsx,
        ngridsy = ngridsy,
        overlap = overlap,
        niterations = niterations,
        px = px,
        py = py,
        tag = tag,
    )
    run_profile(opts)
end

"""
Overlay RSS and CPU time series from one or more `resource_timeseries.csv` files.
Writes a two-panel PNG (RSS and CPU fraction vs elapsed time).
CPU fraction uses each run's `benchmark_results.json` (`metadata.reference_cores`)
when needed; pass `--reference-cores` to override for every series.

# Arguments

- `csv_paths`: paths to `resource_timeseries.csv` files

# Options

- `--labels <list>`: comma-separated labels (default: parent folder names)
- `--reference-cores <n>`: optional override for `cpu_cores_used` normalisation
- `--output <path>`: output PNG path (default: `benchmarks/output/resource_comparison.png`)
"""
@cast function plot(
    csv_paths::String...;
    labels::String = "",
    reference_cores = nothing,
    output::String = "",
)
    labs = isempty(labels) ? String[] : String[strip(s) for s in split(labels, ',') if !isempty(strip(s))]
    out = isempty(output) ? joinpath(BENCHMARK_OUTPUT_DIR, "resource_comparison.png") : output
    ref = reference_cores isa Real ? Float64(reference_cores) : nothing
    plot_resource_timeseries(
        collect(String, csv_paths);
        labels = labs,
        reference_cores = ref,
        output = out,
    )
end

# Initialise Comonicon CLI
Comonicon.@main

end
