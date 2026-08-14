module ismip7_16km_synthetic

using MPI
using WAVI

function driver_grid()
    nx = 381
    ny = 381
    nσ = 26
    x0 = -3048000.0
    y0 = -3048000.0
    dx = 16000.0
    dy = 16000.0

    hm = trues(nx, ny)
    # Add simple boundaries to make the mask meaningful
    hm[[1, end], :] .= false
    hm[:, [1, end]] .= false

    u_iszero = falses(nx+1, ny)
    u_iszero[[1, 2, end-1, end], :] .= true

    v_iszero = falses(nx, ny+1)
    v_iszero[:, [1, 2, end-1, end]] .= true

    # Sigma grid from 0 to 1
    σ = collect(range(0.0, 1.0, length=nσ))

    grid = Grid(
        nx = nx,
        ny = ny,
        nσ = nσ,
        x0 = x0,
        y0 = y0,
        dx = dx,
        dy = dy,
        h_mask = hm,
        u_iszero = u_iszero,
        v_iszero = v_iszero,
        σ = σ,
    )
    return grid
end

function driver_run(; folder = "outputs", grid = driver_grid(), spec = WAVI.BasicSpec())
    nx, ny, nσ = grid.nx, grid.ny, grid.nσ

    t_setup = @elapsed begin
        # Create dummy static arrays mimicking the ISMIP7 setup
    bed = fill(-500.0, nx, ny)

    # Weertman drag coefficient
    weertman_c = fill(1.0e4, nx, ny)
    sliding_law = WeertmanSlidingLaw(drag_coefficient = weertman_c)

    # Simple defaults for params so we don't need forcing files
    params = Params(default_thickness = 1000.0, default_temperature = 260.0, accumulation_rate = 0.3)

    # Solver parameters
    solver_params = SolverParams(maxiter_picard = 5)

    model = Model(
        grid = grid,
        bed_elevation = bed,
        params = params,
        sliding_law = sliding_law,
        spec = spec,
        solver_params = solver_params,
        thermo_dynamics = NoThermoDynamics()
    )

    # Timestepping parameters (just a few steps to benchmark the solve)
    dt = 0.1
    end_time = 0.5  # 5 time steps is enough to establish a benchmark
    chkpt_freq = 10.0 # Don't checkpoint

    timestepping_params = TimesteppingParams(
        niter0 = 0,
        dt = dt,
        end_time = end_time,
        chkpt_freq = chkpt_freq,
        step_thickness = true
    )

    outputs = (
        h = "global_fields.gh.h",
        u = "global_fields.gh.u",
        v = "global_fields.gh.v",
        b = "global_fields.gh.b",
        grfrac = "global_fields.gh.grounded_fraction",
        mpi_rank = "global_fields.gh.mpi_rank",
    )

    # Output parameters (use 'folder' parameter passed by harness)
    output_params = OutputParams(
        outputs = outputs,
        output_freq = 0.5,
        output_path = folder,
        output_format = "jld2",
        zip_format = "nc"
    )

    simulation = Simulation(
        model = model,
        timestepping_params = timestepping_params,
        output_params = output_params,
    )

    end # end setup block

    t_solve = @elapsed begin
        run_simulation!(simulation)
    end

    return (simulation = simulation, setup_time = t_setup, solve_time = t_solve)
end

const driver_plot_vars = ["h", "u", "v"]

end
