# MPISpec with child_architecture = GPU() (run under MPI with 2 ranks).
using Test
using MPI
using WAVI
using WAVI.Architectures: architecture, array_type

try
    using CUDA
catch
end

@testset "MPISpec GPU child_architecture" begin
    MPI.Initialized() || MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)
    @test nprocs == 2

    cuda_ext = Base.get_extension(WAVI, :WAVICUDAExt)
    if cuda_ext === nothing || !(gpu_device() isa GPU)
        @info "Skipping MPI+GPU tests (no functional CUDA device)"
        return
    end

    grid = Grid(nx = 12, ny = 8, nσ = 3, dx = 1.0, dy = 1.0)
    local_spec = ThreadedSpec(ngridsx = 2, ngridsy = 1, overlap = 1, niterations = 1)
    @test_throws ArgumentError MPISpec(
        nprocs, 1, 2, grid;
        pou = true,
        local_spec = local_spec,
        child_architecture = GPU(),
    )

    spec = MPISpec(nprocs, 1, 2, grid; pou = true, niterations = 2, child_architecture = GPU())
    @test architecture(spec) isa GPU
    @test child_architecture(spec) isa GPU

    model = Model(
        grid = grid,
        bed_elevation = -500.0 .* ones(grid.nx, grid.ny),
        params = Params(accumulation_rate = 0.1),
        solver_params = SolverParams(maxiter_picard = 1),
        initial_conditions = InitialConditions(initial_thickness = 100.0 .* ones(grid.nx, grid.ny)),
        spec = spec,
    )
    AT = array_type(architecture(spec))
    @test model.fields.gh.h isa AT
    @test model.fields.gu.u isa AT
    if rank == 0
        @test model.spec.global_fields.gh.h isa Array
    end

    th, rh, bh, lh = WAVI.Specs.get_halos(model.spec)
    h = model.fields.gh.h
    h .= -99.0
    h[(1 + lh):(end - rh), (1 + th):(end - bh)] .= rank + 1.0
    WAVI.Specs.halo_exchange!(model; fields = [:h])
    h_host = Array(h)
    if model.spec.left > -1 && lh > 0
        @test all(h_host[1:lh, :] .== (model.spec.left + 1.0))
    end
    if model.spec.right > -1 && rh > 0
        @test all(h_host[(end - rh + 1):end, :] .== (model.spec.right + 1.0))
    end

    model.fields.gh.h .= -1.0
    model.fields.gh.h[(1 + lh):(end - rh), (1 + th):(end - bh)] .= rank + 10.0
    gathered = WAVI.Specs.collect_mpi_field!(model, [:global_fields, :gh, :h])
    if rank == 0
        vals = Set(vec(gathered))
        @test 10.0 in vals
        @test 11.0 in vals
        @test gathered isa Array
    end

    u_core, v_core = WAVI.Specs.core_inner_masks(model)
    @test u_core isa Vector{Bool}
    @test v_core isa Vector{Bool}
    @test length(u_core) == model.fields.gu.ni
    @test length(v_core) == model.fields.gv.ni

    # Hits the Schwarz residual check (host masks vs GPU packed residual).
    update_velocities!(model)
    @test all(isfinite, Array(model.fields.gu.u))
    @test all(isfinite, Array(model.fields.gv.v))
end

MPI.Finalize()
