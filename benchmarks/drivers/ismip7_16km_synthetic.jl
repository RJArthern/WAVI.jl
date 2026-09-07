module ismip7_16km_synthetic

using MPI
using WAVI

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
const ρ_ice  = 918.0    # kg/m³
const ρ_sw   = 1028.0   # kg/m³

# ---------------------------------------------------------------------------
# MISMIP+-style bed formula (Asay-Davis et al. 2016, scaled to domain)
#
# The bed has two components:
#   B_x(x)   – a 6th-order polynomial in x that creates a retrograde slope
#               (bed deepens from the grounding-line inward, driving marine
#                ice-sheet instability physics)
#   B_y(x,y) – a lateral channel perturbation that carves N_STREAMS ice-stream
#               troughs across the domain width, breaking y-symmetry and
#               coupling u and v throughout the solve
#
# Parameters are adapted from MISMIP+ (640 km × 80 km) to the 6096 km domain
# by preserving the non-dimensional shape and scaling x_bar / wc / fc to match.
# ---------------------------------------------------------------------------

const B0 =  -150.0    # m  (bed elevation offset)
const B2 =  -728.8    # m  (polynomial coefficients from Asay-Davis et al.)
const B4 =   343.91   # m
const B6 =   -50.57   # m
const DC   =   500.0  # m  (channel depth relative to ridge)
const N_STREAMS = 6   # number of ice-stream channels tiled across y

"""
    mismip_bed(x_m, y_m, Lx_m, Ly_m) -> Float64

Return the MISMIP+-style bed elevation (m) at physical coordinates
(x_m, y_m) given domain half-widths Lx_m and Ly_m.
"""
function mismip_bed(x_m::Float64, y_m::Float64, Lx_m::Float64, Ly_m::Float64)
    # Non-dimensionalise x by characteristic length (≈ 50% of domain half-width)
    x_bar = 0.5 * Lx_m
    ξ = x_m / x_bar

    # Polynomial retrograde slope (same non-dimensional shape as MISMIP+)
    Bx = B0 + B2*ξ^2 + B4*ξ^4 + B6*ξ^6

    # Channel half-width and edge width, scaled to domain
    wc = Ly_m / (N_STREAMS * 2.5)   # half-width of each stream channel
    fc = wc / 6.0                    # transition length (sharp channel edges)

    # Tile N_STREAMS evenly spaced channels across full y-extent
    channel_spacing = 2.0 * Ly_m / N_STREAMS
    By = 0.0
    for k in 0:(N_STREAMS - 1)
        y_centre = -Ly_m + (k + 0.5) * channel_spacing
        By += DC / (1.0 + exp(-2.0 * (y_m - y_centre - wc) / fc)) +
              DC / (1.0 + exp( 2.0 * (y_m - y_centre + wc) / fc))
    end
    # Normalise so channels sum to at most DC total depth
    By *= DC / (N_STREAMS * DC)

    return max(Bx + By, -1500.0)  # floor at -1500 m (continental shelf depth)
end

"""
    mismip_drag(x_m, y_m, Lx_m, Ly_m) -> Float64

Return the Weertman drag coefficient (Pa m^{-1/m} s^{1/m}) at (x_m, y_m).
Low drag (10^2) is placed inside each ice-stream channel to match the bed
geometry; high drag (10^4) sits on the inter-stream ridges.
"""
function mismip_drag(x_m::Float64, y_m::Float64, Lx_m::Float64, Ly_m::Float64)
    wc = Ly_m / (N_STREAMS * 2.5)
    fc = wc / 6.0

    channel_spacing = 2.0 * Ly_m / N_STREAMS
    in_channel = 0.0
    for k in 0:(N_STREAMS - 1)
        y_centre = -Ly_m + (k + 0.5) * channel_spacing
        # Smooth logistic mask: 1 inside channel, 0 on ridges
        in_channel += (1.0 / (1.0 + exp( 4.0 * (y_m - y_centre - wc) / fc))) *
                      (1.0 / (1.0 + exp(-4.0 * (y_m - y_centre + wc) / fc)))
    end
    in_channel = clamp(in_channel, 0.0, 1.0)

    C_ridge  = 1.0e4   # sticky inter-stream ridge
    C_stream = 1.0e2   # slippery ice-stream channel
    return C_ridge * (1.0 - in_channel) + C_stream * in_channel
end

# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------

function driver_grid()
    nx = 381
    ny = 381
    nσ = 26
    x0 = -3048000.0   # domain centred at origin: x ∈ [-3048, +3048] km
    y0 = -3048000.0   # y ∈ [-3048, +3048] km
    dx = 16000.0
    dy = 16000.0

    hm = trues(nx, ny)
    hm[[1, end], :] .= false
    hm[:, [1, end]] .= false

    u_iszero = falses(nx+1, ny)
    u_iszero[[1, 2, end-1, end], :] .= true

    v_iszero = falses(nx, ny+1)
    v_iszero[:, [1, 2, end-1, end]] .= true

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
    nx, ny = grid.nx, grid.ny

    t_setup = @elapsed begin
        dx, dy = grid.dx, grid.dy
        x0, y0 = grid.x0, grid.y0

        # Physical cell-centre coordinates (m)
        xs = [x0 + (i - 0.5) * dx for i in 1:nx]
        ys = [y0 + (j - 0.5) * dy for j in 1:ny]
        Lx_m = nx * dx / 2.0
        Ly_m = ny * dy / 2.0

        # ------------------------------------------------------------------
        # Bed elevation – MISMIP+-style retrograde slope + lateral channels
        # ------------------------------------------------------------------
        bed = [mismip_bed(xs[i], ys[j], Lx_m, Ly_m) for i in 1:nx, j in 1:ny]

        # ------------------------------------------------------------------
        # Basal drag – low inside ice-stream channels, high on ridges
        # ------------------------------------------------------------------
        weertman_c = [mismip_drag(xs[i], ys[j], Lx_m, Ly_m) for i in 1:nx, j in 1:ny]
        sliding_law = WeertmanSlidingLaw(drag_coefficient = weertman_c)

        # ------------------------------------------------------------------
        # Initial thickness – parabolic dome to drive outward ice flow.
        # An ice sheet on a continent is thick at the interior divide and
        # thins towards its marine margins.
        # ------------------------------------------------------------------
        H_max = 3000.0
        H_min =  500.0
        initial_thickness = [
            H_min + (H_max - H_min) * max(0.0, 1.0 - (xs[i]^2 + ys[j]^2) / (0.7 * Lx_m)^2)
            for i in 1:nx, j in 1:ny
        ]

        # ------------------------------------------------------------------
        # Params and solver
        # ------------------------------------------------------------------
        params = Params(
            default_thickness   = H_max,
            default_temperature = 260.0,
            accumulation_rate   = 0.3,
        )

        # Force exactly 5 Picard iterations to avoid early-convergence skew
        solver_params = SolverParams(maxiter_picard = 5, tol_picard = 0.0)

        model = Model(
            grid               = grid,
            bed_elevation      = bed,
            params             = params,
            sliding_law        = sliding_law,
            spec               = spec,
            solver_params      = solver_params,
            initial_conditions = InitialConditions(initial_thickness = initial_thickness),
            thermo_dynamics    = NoThermoDynamics(),
        )

        # Timestepping: 5 steps, no checkpointing
        dt         = 0.1
        end_time   = 0.5
        chkpt_freq = 10.0

        timestepping_params = TimesteppingParams(
            niter0         = 0,
            dt             = dt,
            end_time       = end_time,
            chkpt_freq     = chkpt_freq,
            step_thickness = true,
        )

        outputs = (
            h        = "global_fields.gh.h",
            u        = "global_fields.gh.u",
            v        = "global_fields.gh.v",
            b        = "global_fields.gh.b",
            grfrac   = "global_fields.gh.grounded_fraction",
            mpi_rank = "global_fields.gh.mpi_rank",
        )

        output_params = OutputParams(
            outputs       = outputs,
            output_freq   = 0.5,
            output_path   = folder,
            output_format = "jld2",
            zip_format    = "nc",
        )

        simulation = Simulation(
            model               = model,
            timestepping_params = timestepping_params,
            output_params       = output_params,
        )
    end # end setup block

    t_solve = @elapsed begin
        run_simulation!(simulation)
    end

    return (simulation = simulation, setup_time = t_setup, solve_time = t_solve)
end

const driver_plot_vars = ["h", "u", "v"]

end
