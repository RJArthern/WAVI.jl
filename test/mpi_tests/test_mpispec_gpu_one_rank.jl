# One-rank MPISpec GPU: halo exchange must skip packing (run under MPI with 1 rank).
using Test
using MPI
using WAVI
using WAVI.Architectures: architecture, array_type

try
    using CUDA
catch
end

@testset "MPISpec one rank skips halo packing" begin
    MPI.Initialized() || MPI.Init()
    comm = MPI.COMM_WORLD
    @test MPI.Comm_size(comm) == 1

    cuda_ext = Base.get_extension(WAVI, :WAVICUDAExt)
    if cuda_ext === nothing || !(gpu_device() isa GPU)
        @info "Skipping one-rank MPI+GPU tests (no functional CUDA device)"
        return
    end

    grid = Grid(nx = 12, ny = 8, nσ = 3, dx = 1.0e3, dy = 1.0e3)
    spec = MPISpec(1, 1, 2, grid; pou = true, niterations = 2, child_architecture = GPU())
    @test !WAVI.Specs.mpi_has_neighbours(spec)

    model = Model(
        grid = grid,
        bed_elevation = -500.0 .* ones(grid.nx, grid.ny),
        params = Params(accumulation_rate = 0.1),
        solver_params = SolverParams(maxiter_picard = 1),
        initial_conditions = InitialConditions(initial_thickness = 100.0 .* ones(grid.nx, grid.ny)),
        spec = spec,
        verbose = false,
    )
    @test model.fields.gh.h isa array_type(architecture(spec))

    WAVI.Specs.halo_exchange!(model; fields = [:h])
    @test model.spec.halo_scratch === nothing
end

MPI.Finalize()
