# Thin Comonicon entry point for the WAVIBenchmarks harness.
#
# From the repository root dir:
#
#   julia --project=benchmarks benchmarks/run.jl drivers
#
#   julia --project=benchmarks -t 1 benchmarks/run.jl run basic mismip_plus
#
#   julia --project=benchmarks -t 4 benchmarks/run.jl run threaded mismip_plus \
#     --ngridsx 2 --ngridsy 2
#
#   mpiexecjl --project=benchmarks -n 4 julia -t 1 benchmarks/run.jl run mpi mismip_plus \
#     --px 4 --py 1
#
#   mpiexecjl --project=benchmarks -n 4 julia -t 1 benchmarks/run.jl run mpi_gpu ismip7_16km_synthetic \
#     --px 4 --py 1
#   (Needs CUDA.jl, a GPU node, and one GPU per rank. Do not run `module load cuda`.)
#
#   (Avoid --px 2 --py 2 on MISMIP+: the y-patches are too thin for efficient PoU.)

# Force GR plotting backend to run headlessly to prevent flashing plot window
ENV["GKSwstype"] = "100"

# Load CUDA at top-level WAVI, else, getting `CUDA.functional()` (and GPUSpec) is too new for the current world age warnings.
if "gpu" in ARGS || "mpi_gpu" in ARGS
    using CUDA
end

using WAVIBenchmarks

if abspath(PROGRAM_FILE) == @__FILE__
    # Reconstruct command used to launch benchmark for reproducibility logging
    # Start with the Julia binary and the active thread count
    parts = String["julia", "-t $(Threads.nthreads())"]
    
    # Include active project path if set
    proj = Base.active_project()
    if !isnothing(proj)
        push!(parts, "--project=$(relpath(dirname(proj), pwd()))")
    end
    
    # Append script name (e.g. 'benchmarks/run.jl') and all CLI arguments
    push!(parts, PROGRAM_FILE)
    append!(parts, ARGS)
    
    # Store reconstructed command string into global reference inside the harness
    WAVIBenchmarks.set_benchmark_command!(join(parts, " "))
    
    # Parse CLI args and delegate subcommand execution to Comonicon
    exit(WAVIBenchmarks.command_main(ARGS))
end
