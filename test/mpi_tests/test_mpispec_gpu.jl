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

    cpu_spec = MPISpec(nprocs, 1, 2, grid; pou = true, niterations = 2)
    cpu_model = Model(
        grid = grid,
        bed_elevation = -500.0 .* ones(grid.nx, grid.ny),
        params = Params(accumulation_rate = 0.1),
        solver_params = SolverParams(maxiter_picard = 1),
        initial_conditions = InitialConditions(initial_thickness = 100.0 .* ones(grid.nx, grid.ny)),
        spec = cpu_spec,
    )
    cpu_h = cpu_model.fields.gh.h
    cpu_h .= -99.0
    cpu_h[(1 + lh):(end - rh), (1 + th):(end - bh)] .= rank + 1.0
    WAVI.Specs.halo_exchange!(cpu_model; fields = [:h])
    @test h_host ≈ cpu_h rtol = 1e-8 atol = 1e-8

    scratch = WAVI.Specs.ensure_mpi_halo_scratch!(model)
    @test scratch.dev_halo_strip !== nothing
    @test !(scratch.dev_halo_strip isa Array)
    field = AT(reshape(collect(1.0:300.0), 20, 15))
    buf = Float64[]
    packed_x = WAVI.Specs.pack_halo_strip!(buf, field, 2:3, 1:15, scratch)
    @test reshape(collect(packed_x), 2, 15) ≈ Array(field)[2:3, :]
    packed_y = WAVI.Specs.pack_halo_strip!(buf, field, 1:20, 4:5, scratch)
    @test reshape(collect(packed_y), 20, 2) ≈ Array(field)[:, 4:5]

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
    @test u_core isa AbstractVector{Bool}
    @test v_core isa AbstractVector{Bool}
    @test !(u_core isa Vector{Bool})
    @test !(v_core isa Vector{Bool})
    gu, gv = model.fields.gu, model.fields.gv
    @test length(u_core) == gu.ni
    @test length(v_core) == gv.ni

    # Non-zero packed vector: device reduce must match the host core-only formula.
    n = gu.ni + gv.ni
    packed = similar(gu.u, Float64, n)
    copyto!(packed, Float64.(1:n))
    u_range = 1:gu.ni
    v_range = (gu.ni + 1):n
    gpu_sq = WAVI.Specs.masked_sum_abs2(packed, u_range, u_core) +
             WAVI.Specs.masked_sum_abs2(packed, v_range, v_core)
    packed_h = Array(packed)
    u_m, v_m = Array(u_core), Array(v_core)
    host_sq = sum(abs2, packed_h[u_range][u_m]) + sum(abs2, packed_h[v_range][v_m])
    @test host_sq > 0
    @test gpu_sq ≈ host_sq rtol = 1e-8 atol = 1e-8

    # Same ice as the CPU halo check; collect above overwrote GPU thickness.
    copyto!(model.fields.gh.h, cpu_h)
    WAVI.Processes.inner_update!(cpu_model)
    WAVI.Processes.inner_update!(model)

    function core_sq(m)
        s = WAVI.Utilities.stencil_scratch!(m)
        resid = s.picard_resid
        rhs = WAVI.get_rhs(m)
        WAVI.get_resid!(resid, WAVI.get_start_guess(m), WAVI.get_op(m), rhs)
        uc, vc = WAVI.Specs.core_inner_masks(m)
        ni, nv = m.fields.gu.ni, m.fields.gv.ni
        rsq = WAVI.Specs.masked_sum_abs2(resid, 1:ni, uc) +
              WAVI.Specs.masked_sum_abs2(resid, (ni + 1):(ni + nv), vc)
        bsq = WAVI.Specs.masked_sum_abs2(rhs, 1:ni, uc) +
              WAVI.Specs.masked_sum_abs2(rhs, (ni + 1):(ni + nv), vc)
        return rsq, bsq
    end
    gpu_resid_sq, gpu_rhs_sq = core_sq(model)
    cpu_resid_sq, cpu_rhs_sq = core_sq(cpu_model)
    @test gpu_resid_sq ≈ cpu_resid_sq rtol = 1e-8 atol = 1e-8
    @test gpu_rhs_sq ≈ cpu_rhs_sq rtol = 1e-8 atol = 1e-8

    # Hits the Schwarz residual check (device reduce, Allreduce two scalars).
    update_velocities!(model)
    @test all(isfinite, Array(model.fields.gu.u))
    @test all(isfinite, Array(model.fields.gv.v))

    @testset "MPISpec GPU checkpoint pickup and resume" begin
        dt = 0.1
        t_mid = 0.2
        t_end = 0.4
        mktempdir() do dir
            function gpu_mpi_model()
                spec = MPISpec(nprocs, 1, 2, grid; pou = true, niterations = 2, child_architecture = GPU())
                return Model(
                    grid = grid,
                    bed_elevation = -500.0 .* ones(grid.nx, grid.ny),
                    params = Params(accumulation_rate = 0.1),
                    solver_params = SolverParams(maxiter_picard = 1),
                    initial_conditions = InitialConditions(initial_thickness = 100.0 .* ones(grid.nx, grid.ny)),
                    spec = spec,
                    verbose = false,
                )
            end
            sim_write = Simulation(
                model = gpu_mpi_model(),
                timestepping_params = TimesteppingParams(
                    dt = dt,
                    end_time = t_mid,
                    chkpt_freq = t_mid,
                    chkpt_path = dir,
                ),
                output_params = OutputParams(output_path = dir),
            )
            run_simulation!(sim_write)
            n_iter = sim_write.clock.n_iter
            @test isfile(joinpath(dir, WAVI.Outputs.checkpoint_filename(sim_write.model.spec, n_iter)))
            h0 = Array(sim_write.model.fields.gh.h)
            u0 = Array(sim_write.model.fields.gu.u)

            sim_pick = Simulation(
                model = gpu_mpi_model(),
                timestepping_params = TimesteppingParams(
                    dt = dt,
                    end_time = t_end,
                    niter0 = n_iter,
                    chkpt_path = dir,
                ),
                output_params = OutputParams(output_path = dir),
            )
            @test sim_pick.clock.n_iter == n_iter
            @test sim_pick.model.fields.gh.h isa AT
            @test sim_pick.model.fields.stencil_scratch[] === nothing
            @test Array(sim_pick.model.fields.gh.h) ≈ h0
            @test Array(sim_pick.model.fields.gu.u) ≈ u0

            run_simulation!(sim_pick)
            @test sim_pick.clock.time ≈ t_end
            @test sim_pick.model.fields.gh.h isa AT
            @test all(isfinite, Array(sim_pick.model.fields.gu.u))

            sim_ctrl = Simulation(
                model = gpu_mpi_model(),
                timestepping_params = TimesteppingParams(
                    dt = dt,
                    end_time = t_end,
                    chkpt_path = dir,
                ),
                output_params = OutputParams(output_path = dir),
            )
            run_simulation!(sim_ctrl)
            @test Array(sim_pick.model.fields.gh.h) ≈ Array(sim_ctrl.model.fields.gh.h)
            @test Array(sim_pick.model.fields.gu.u) ≈ Array(sim_ctrl.model.fields.gu.u)
        end
    end
end

MPI.Finalize()
