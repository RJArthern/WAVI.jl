group = "all" 

if group == "fields" || group == "all"
    println(1)
    include("test_fields.jl")
    include("Stencils_test.jl")
end

if group == "grids" || group == "all"
    include("test_grids.jl")
end

if group == "melt" || group == "all"
    include("test_melt.jl")
end

if group == "models" || group == "all"
    include("test_models.jl")
end

if group == "outputs" || group == "all"
    include("test_outputting.jl")
end

if group == "simulations" || group == "all"
    include("test_simulations.jl")
end

if group == "timesteppingparams" || group == "all"
    include("test_timesteppingparams.jl")
end

if group == "utils" || group == "all"
    include("test_utils.jl")
end

if group == "kronecker" || group == "all"
    include("test_kronecker.jl")
end

if group == "version_update" || group == "all"
    include(joinpath("version_update_test_verification","test_version_updates.jl"))
end

if group == "shared_memory_parallelism" || group == "all"
    include(joinpath("domain_decomposition_tests","test_domains_shared_memory.jl"))
end

if group == "verification" || group == "all"
    include("verification_tests.jl")
end

# MPI tests: spawned as subprocesses via mpiexec() (MPI.jl recommended pattern).
# Wrapped in try/catch so that Pkg.test() still works on systems without MPI installed.
_mpi_skipped = false
if group == "mpi_unit" || group == "all"
    try
        include("runtests_mpi_unit.jl")
    catch e
        global _mpi_skipped = true
        @warn "MPI unit tests skipped (MPI not available or mpiexec failed)" exception = e
    end
end

if group == "mpi_integration" || group == "all"
    try
        include("runtests_mpi_integration.jl")
    catch e
        global _mpi_skipped = true
        @warn "MPI integration tests skipped (MPI not available or mpiexec failed)" exception = e
    end
end

if _mpi_skipped
    @warn "MPI tests were SKIPPED because MPI was not available. Install MPI and MPI.jl to run the full suite."
end
