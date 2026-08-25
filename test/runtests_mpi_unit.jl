using MPI
using Test
using WAVI

try
    using CUDA
catch
end

include("run_mpi_script.jl")

@testset "MPISpec MPI unit tests" begin
    p, cmd = run_mpi_script("test_mpispec_halo_collect.jl"; nprocs = 2)
    @test success(p)
    if !success(p)
        @error "MPI unit tests failed" cmd exitcode = p.exitcode
    end
end

@testset "MPISpec GPU child_architecture" begin
    if Base.get_extension(WAVI, :WAVICUDAExt) === nothing || !(gpu_device() isa GPU)
        @info "Skipping MPI+GPU unit tests (no functional CUDA device)"
        @test_skip "MPI+GPU (no functional CUDA device)"
    else
        p, cmd = run_mpi_script("test_mpispec_gpu.jl"; nprocs = 2)
        @test success(p)
        if !success(p)
            @error "MPI+GPU unit tests failed" cmd exitcode = p.exitcode
        end
    end
end
