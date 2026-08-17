export MPISpec

using JLD2
using MPI
using Parameters

using WAVI.Parameters

import WAVI: AbstractField, AbstractGrid, AbstractMeltRate, AbstractSurfaceMassBalance, AbstractFracture, AbstractSlidingLaw, 
                   AbstractBasalHydrology, AbstractThermoDynamics, AbstractModel
import WAVI.Deferred: Collector, register_item!, field_extractor
import WAVI.Fields: GridField, InitialConditions, HGrid, UGrid, VGrid, CGrid, SigmaGrid
import WAVI.Grids: Grid
import WAVI.MeltRates: UniformMeltRate
import WAVI.Models: BasicSpec, Model, get_bed_elevation
import WAVI.Outputs: write_outputs, zip_output, OutputParams, checkpoint_filename, load_checkpoint, write_checkpoint!, checkpoint_path, with_cleared_stencil_scratch
import WAVI.Parameters: TimesteppingParams
import WAVI.Processes: update_state!, update_model_velocities!, update_velocities!, update_velocities_on_h_grid!, inner_update!, precondition!, update_preconditioner!, update_rheological_operators!
import WAVI.Simulations: run_simulation!, timestep!, update_model_climate_forcing!
import WAVI.Time: Clock
import WAVI.Wavelets: UWavelets, VWavelets

# FIXME: important to realise that this specification has become a complex structure to house many things that should be baked into the model structurally
#  not least the global grid and fields 

"""
Reusable PoU weights, velocity workspaces, and strip packs for neighbour prolong.
Allocated once per local velocity shape; reused across Schwarz iterations.
"""
mutable struct MPIPoUScratch{T <: AbstractFloat}
    ωu::Matrix{T}
    ωv::Matrix{T}
    work_u::Matrix{T}
    work_v::Matrix{T}
    send_l::Vector{T}
    recv_l::Vector{T}
    send_r::Vector{T}
    recv_r::Vector{T}
    send_t::Vector{T}
    recv_t::Vector{T}
    send_b::Vector{T}
    recv_b::Vector{T}
end

function MPIPoUScratch(ωu::Matrix{T}, ωv::Matrix{T}) where {T <: AbstractFloat}
    empty = T[]
    return MPIPoUScratch{T}(
        ωu, ωv, similar(ωu), similar(ωv),
        copy(empty), copy(empty), copy(empty), copy(empty),
        copy(empty), copy(empty), copy(empty), copy(empty),
    )
end

"""
Struct to represent the MPI parallel specification of a model.

Fields:
    px, py: Number of processes in x and y directions
    halo: Number of halo cells on each side of the subgrid
    rank: MPI rank of the current process
    comm: MPI communicator for the current process
    coords: Cartesian coordinates of the current process
    top, right, bottom, left: Neighbouring processes in the x and y directions
    global_size: Total number of processes
    global_comm: MPI communicator for the global grid
    global_grid: Global grid to be used for the model
    damping: Damping factor for halo exchange (default=0.0)
    niterations: Number of Schwarz iterations per Picard iteration (default=5)
    field_collector: Field collector for the model
    pou_scratch: Cached PoU weights / strip buffers (filled on first prolong)
    local_spec: Optional intraprocess ThreadedSpec for the rank-local solve (default=nothing to wavelet)
"""
mutable struct MPISpec{N <: Integer, T <: Number, M <: MPI.Comm, G <: AbstractGrid} <: AbstractDecompSpec
    # MPI Specification information
    px::N
    py::N
    halo::N

    # MPI process information
    global_size::N
    global_comm::M
    global_grid::G
    global_fields::Union{GridField, Nothing}
    field_collector::Collector

    rank::N
    comm::M
    coords::Array{N, 1}

    # Neighbourhood information
    top::Union{N, Nothing}
    right::Union{N, Nothing}
    bottom::Union{N, Nothing}
    left::Union{N, Nothing}
    pou::Bool
    damping::T
    niterations::N  # Number of Schwarz iterations per Picard iteration
    pou_scratch::Union{Nothing, MPIPoUScratch}
    local_spec::Union{Nothing, ThreadedSpec}

    @doc """
    Constructor for MPISpec.

    px, py: Number of processes in x and y directions
    halo: Number of halo cells on each side of the subgrid
    grid: Grid to be used for the model
    niterations: Number of Schwarz iterations per Picard iteration (default=5)
    pou: Whether to use partition of unity for halo exchange (default=true)
    damping: Damping factor for halo exchange (default=0.0)
    local_spec: Optional ThreadedSpec for the intraprocess local solve (default=nothing)
    """
    function MPISpec(px::Integer, py::Integer, halo::Integer, grid::AbstractGrid; pou::Bool=true, damping::AbstractFloat=0.0, niterations::Integer=5, local_spec::Union{Nothing,ThreadedSpec}=nothing)
        (px < 1 || py < 1 || halo < 0) && 
            throw(ArgumentError("Invalid parameters specified for MPISpec"))

        MPI.Initialized() || MPI.Init()
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        size = MPI.Comm_size(comm)
        @debug "Creating dimensions of $(size) with ($(px), $(py))"

        if size > 1 && halo == 0
            throw(ArgumentError(
                "MPI Configuration Error: halo must be >= 1 when running with multiple MPI ranks (size=$(size))."
            ))
        end

        if size > 1 && pou && halo < 2
            throw(ArgumentError(
                "MPI Configuration Error: Partition of Unity (pou=true) requires a minimum halo size of 2 for blending."
            ))
        end

        # Create Virtual Cartesian Topology based on no. of procs in each direction
        dims = MPI.Dims_create(size, (px, py))
        cart_comm = MPI.Cart_create(comm, dims)
        x_coord, y_coord = MPI.Cart_coords(cart_comm)

        # Return the rank of the process in the direction specified
        # Note: MPI_Cart_shift returns `MPI_PROC_NULL` (represented as a negative integer)
        #       if the neighbour is not present
        top = MPI.Cart_shift(cart_comm, 1, -1)[2]
        right = MPI.Cart_shift(cart_comm, 0, 1)[2]
        bottom = MPI.Cart_shift(cart_comm, 1, 1)[2]
        left = MPI.Cart_shift(cart_comm, 0, -1)[2]

        @debug "[$(rank+1)/$(size)] Neighbours $(top),$(right),$(bottom),$(left)"

        # Check if the process grid is valid
        if px * py != size
            valid = [(i, div(size,i)) for i in 1:size if size % i == 0]
            error(
                "MPI Configuration Error: Invalid process grid.\n\n" *
                "size = $(size), px = $(px), py = $(py)\n\n" *
                "px * py must equal size.\n\n" *
                "Valid (px, py) combinations for size=$(size):\n" *
                join(["  ($(a), $(b))" for (a,b) in valid], "\n")
            )
        end

        # Validation checks for grid size vs halo
        # Need to ensure local grid size is sufficient for halo exchange
        # The constraint is generally: `Core > Halo` (for edge nodes) or `Core > 2*Halo` (for interior nodes)
        # where `Core = div(global_dim, procs)`
        if size != 1
            validate_dimension("nx", grid.nx, px, halo)
            validate_dimension("ny", grid.ny, py, halo)
        end

        return new{Int, Float64, MPI.Comm, AbstractGrid}(
            px, 
            py, 
            halo,
            size,
            comm,
            grid,
            nothing,
            Collector(),
            rank,
            cart_comm,
            [x_coord, y_coord],
            top, right, bottom, left,
            pou,
            damping,
            niterations,
            nothing,  # pou_scratch filled on first PoU prolong
            local_spec,
            )
    end
end

include("MPI/utils.jl")
include("MPI/exchanges.jl")
include("MPI/outputs.jl")
include("MPI/mpi_checkpoints.jl")

function Model(grid::G,
               bed_elevation::Union{Integer, Function, AbstractArray},
               spec::S;
               initial_conditions::InitialConditions = InitialConditions(),
               params::Params = Params(),
               solver_params::SolverParams = SolverParams(),
               shelf_melt_rate::M = UniformMeltRate(),
               surface_mass_balance::SMB = AccumulationFromParams(),
               fracture::FR = ConstantDamage(),
               sliding_law::SL = WeertmanSlidingLaw(),
               basal_hydrology::BH = NoHydrology(),
               thermo_dynamics::TD = NoThermoDynamics(),
               verbose = true)                   where {G<:AbstractGrid, 
                                                        S<:MPISpec, 
                                                        M<:AbstractMeltRate,
                                                        SMB<:AbstractSurfaceMassBalance,
                                                        FR<:AbstractFracture,
                                                        SL<:AbstractSlidingLaw,
                                                        BH<:AbstractBasalHydrology,
                                                        TD<:AbstractThermoDynamics}

    @unpack coords, global_size, rank, comm = spec

    # Recalculate grid dimensions and mask parameters, creating a new local Grid
    th, rh, bh, lh = get_halos(spec)
    nx_local, ny_local = get_size(spec)
    x_start, x_end, y_start, y_end = get_bounds(spec)
    x0_local = grid.x0 + (x_start-1) * grid.dx
    y0_local = grid.y0 + (y_start-1) * grid.dy
    @info "[$(rank+1)/$(global_size)] proc $(coords[1]),$(coords[2]) grid $(nx_local)x$(ny_local) X[$(x_start):$(x_end)] Y[$(y_start):$(y_end)]"
    @debug "[$(rank+1)/$(global_size)] - $(coords) - creating Grid and Model for MPI rank $(rank)"
    @debug "[$(rank+1)/$(global_size)] - proc $(coords[1]),$(coords[2]) - $(th), $(rh), $(bh), $(lh)"
    @debug "[$(rank+1)/$(global_size)] - proc $(coords[1]),$(coords[2]) - grid $(nx_local)x$(ny_local)"
    @debug "[$(rank+1)/$(global_size)] - X [$(x_start):$(x_end)] - Y [$(y_start):$(y_end)] - Centroid $(x0_local),$(y0_local) "

    u_grid_size, v_grid_size = (grid.nx+1, grid.ny), (grid.nx, grid.ny+1)
    
    #expand scalar paramaters onto grid
    params = reconstruct_on_grid(params,grid)

    #Replace all NaN entries with defaults from params on correct grid              
    initial_conditions = reconstruct_on_grid(initial_conditions, params, grid)

    #expand spatial parameters onto grid
    shelf_melt_rate = reconstruct_on_grid(shelf_melt_rate,grid)
    surface_mass_balance = reconstruct_on_grid(surface_mass_balance,grid)
    fracture = reconstruct_on_grid(fracture,grid)
    sliding_law = reconstruct_on_grid(sliding_law,grid)
    basal_hydrology = reconstruct_on_grid(basal_hydrology,grid)
    thermo_dynamics = reconstruct_on_grid(thermo_dynamics,grid)

    #trim initial conditions to local domain
    local_initial_conditions = reconstruct_on_subdomain(initial_conditions, grid, (x_start, x_end, y_start, y_end))

    # dt cannot be copied via the external constructor so we create the structure directly
    local_params = reconstruct_on_subdomain(params,grid,(x_start, x_end, y_start, y_end))

    local_shelf_melt_rate = reconstruct_on_subdomain(shelf_melt_rate,grid,(x_start, x_end, y_start, y_end))
    local_surface_mass_balance = reconstruct_on_subdomain(surface_mass_balance,grid,(x_start, x_end, y_start, y_end))
    local_fracture = reconstruct_on_subdomain(fracture,grid,(x_start, x_end, y_start, y_end))
    local_sliding_law = reconstruct_on_subdomain(sliding_law,grid,(x_start, x_end, y_start, y_end))
    local_basal_hydrology = reconstruct_on_subdomain(basal_hydrology,grid,(x_start, x_end, y_start, y_end))
    local_thermo_dynamics = reconstruct_on_subdomain(thermo_dynamics,grid,(x_start, x_end, y_start, y_end))

    
    u_isfixed = grid.u_isfixed[x_start:x_end+1, y_start:y_end]
    v_isfixed = grid.v_isfixed[x_start:x_end, y_start:y_end+1]
    # RAS/Schwarz: Only fix the outermost edge cell - this provides Dirichlet BC from neighbour
    # Interior halo cells (2:halo) are SOLVED locally, then DISCARDED during exchange
    # The simple overwrite in halo_exchange! replaces them with neighbour's core values
    (spec.left < 0) || (u_isfixed[1,:] .= true; v_isfixed[1,:] .= true)
    (spec.right < 0) || (u_isfixed[end,:] .= true; v_isfixed[end,:] .= true)
    (spec.top < 0) || (u_isfixed[:,1] .= true; v_isfixed[:,1] .= true)
    (spec.bottom < 0) || (u_isfixed[:,end] .= true; v_isfixed[:,end] .= true)
    
    #TODO Halo for hyd_potential_isfixed

    local_grid = Grid(
        nx = nx_local,
        ny = ny_local,
        nσ = grid.nσ,
        dx = grid.dx,
        dy = grid.dy,
        x0 = x0_local,
        y0 = y0_local,
        h_mask = grid.h_mask[x_start:x_end, y_start:y_end],
        h_isfixed = grid.h_isfixed[x_start:x_end, y_start:y_end],
        hyd_potential_isfixed = grid.hyd_potential_isfixed[x_start:x_end, y_start:y_end],
        u_iszero = grid.u_iszero[x_start:x_end+1, y_start:y_end],
        v_iszero = grid.v_iszero[x_start:x_end, y_start:y_end+1],
        u_isfixed = u_isfixed,
        v_isfixed = v_isfixed,
        quadrature_weights = grid.quadrature_weights,
        σ = grid.σ,
        basin_ID = grid.basin_ID[x_start:x_end, y_start:y_end])

    if typeof(bed_elevation) <: AbstractArray
        bed_array = bed_elevation[x_start:x_end, y_start:y_end]
    else
        bed_array = get_bed_elevation(bed_elevation, local_grid)
    end

    # Create local mpi_rank field filled with this rank's number
    local_mpi_rank = fill(Float64(rank), nx_local, ny_local)

    fields = GridField(local_grid, bed_array; initial_conditions=local_initial_conditions, params=local_params, solver_params, mpi_rank=local_mpi_rank)
    model = Model(local_grid, fields, local_params, solver_params, spec, local_shelf_melt_rate, local_surface_mass_balance, local_fracture, local_sliding_law, local_basal_hydrology, 
    local_thermo_dynamics, verbose)

    global_bed = typeof(bed_elevation) <: AbstractArray ? bed_elevation : get_bed_elevation(bed_elevation, grid)
    # Create global mpi_rank field (will be populated during collection)
    global_mpi_rank = zeros(Float64, grid.nx, grid.ny)

    # Global assembly buffers for gather/Bcast only: sized arrays, no IC/operator cost.
    # (Physics lives on each rank's local GridField; bed is retained for static outputs.)
    model.spec.global_fields = GridField(grid, global_bed; initial_conditions, params, solver_params, mpi_rank=global_mpi_rank, assembly_buffer=true)

    MPI.Barrier(comm)
    if rank == 0
        @info "MPI: local models ready on $(global_size) ranks"
    end
    return model
end

##
# Make Model interface work transparently with MPI distributed elements
#

# TODO: remove getproperty, if it's not a global registered field it should be ignored
function Base.getproperty(model::Model{T,N,<:MPISpec,F,G,M,SMB,FR,SL,BH,TD}, s::Symbol) where {T,N,F,G,M,SMB,FR,SL,BH,TD}
    if s == :global_fields
        # TODO: these need to be registered fields, not user-specified
        ## TODO: fields = collate_global_fields(model.fields, model.spec)
        return model.spec.global_fields
    elseif s == :global_grid
        return model.spec.global_grid
    end
    return getfield(model, s)
end

##
# Override Model oriented methods to intercept calls that need extra processing
#

function update_model_velocities!(model::Model{<:Any, <:Any, <:MPISpec})
    update_velocities!(model)
    return model
end

function update_preconditioner!(model::Model{<:Any, <:Any, <:MPISpec})
    ls = model.spec.local_spec
    ls !== nothing && update_preconditioner!(model, ls)
    return model
end

# Rank-local solve: ThreadedSpec when nested, otherwise the default wavelet path.
function local_precondition!(model::Model{<:Any, <:Any, <:MPISpec})
    ls = model.spec.local_spec
    if ls !== nothing
        precondition!(model, ls)
    else
        invoke(precondition!, Tuple{AbstractModel}, model)
    end
end

##
# Overrides for other functionalities that need restricting to root node
#
# TODO: override @debug, @info, @warn and @error for MPI based logging, with the rank out of size and / or grid location

##
# Implementations that affect the simuation and data collection
#
function timestep!(model::AbstractModel{T,N,S},
                   timestepping_params::TimesteppingParams,
                   output_params::OutputParams,
                   clock::Clock) where {T,N,S<:MPISpec}
    if mod(clock.n_iter, timestepping_params.ntimesteps_climate_forcing_update) == 0
        update_model_climate_forcing!(model, clock)
    end

    update_state!(model, clock)

    #write solution if at the first timestep (hack for https://github.com/RJArthern/WAVI.jl/issues/46 until synchronicity is fixed)
    # Have made the interface consistent
    # Have also removed the dependence on individual call
    if (output_params.output_start) && (clock.n_iter == 0)
        collect!(model.spec.field_collector, model)
        write_outputs(model, timestepping_params, output_params, clock)
        clear!(model.spec.field_collector)
    end
    
    if timestepping_params.step_thickness
        update_thickness!(model, timestepping_params)
        # Sync thickness halos; do not RAS-overwrite u/v when using PoU (see mpi_sync_halos_after_thickness!)
        mpi_sync_halos_after_thickness!(model)
    end
    update_clock!(clock, timestepping_params)

    # Collect AFTER thickness update to match BasicSpec output timing
    collect!(model.spec.field_collector, model)
    write_outputs(model, timestepping_params, output_params, clock)
    clear!(model.spec.field_collector)
end

function run_simulation!(model::AbstractModel{T,N,S}, 
                         timestepping_params::TimesteppingParams, 
                         output_params::OutputParams,
                         clock::Clock) where {T,N,S<:MPISpec}
    for field in values(output_params.outputs.items)
        if model.spec.rank == 0
            @debug "Registering $(field.path) from outputs"
        end
        if field.path[1] == :global_fields
            register_mpi_field!(model.spec.field_collector, field.path)
        end
    end

    # TODO: we potentially register other fields here too, but currently concentrating on outputs (update_thickness might want to exploit this mechanism)

    rank = model.spec.rank
    if rank == 0
        @info "MPI: exchanging initial halos"
    end
    mpi_sync_halos_initial!(model)

    for i = (clock.n_iter+1):timestepping_params.n_iter_total
        if rank == 0
            @info "Running iteration $(clock.n_iter)/$(timestepping_params.n_iter_total)"
            if clock.n_iter == 0
                @info "MPI: first velocity solve compiles kernels and may take a while"
            end
        end
        timestep!(model, timestepping_params, output_params, clock)
    end

    zip_output(model, output_params)
end

function register_mpi_field!(collector::Collector, path::Vector{Symbol})
    accessor = function(model)
        collect_mpi_field!(model, path)
        result = model.spec
        for field in path
            result = getproperty(result, field)
        end
        return result
    end
    extractor = field_extractor(join(string.(path), "."), accessor, path)
    register_item!(collector, extractor)
end

"""
    inner_update!(model::Model{<:Any, <:Any, <:MPISpec})

Overload for the inner update of the velocity solve.
We need to sync halos before performing the update so that the viscosity and other
calculations have correct boundary information from neighbouring procs.

Note: `update_rheological_operators!` is called a second time here (after the halo sync)
because the base `inner_update!` builds the operators before we have correct cross-rank
rheology values (β, βeff, ηav). The second call rebuilds them with consistent halo data.
"""
function inner_update!(model::Model{<:Any, <:Any, <:MPISpec})
    # RAS ghost sync for velocities; AS-PoU keeps overlap values from the last prolongation.
    if !model.spec.pou
        halo_exchange!(model; fields=[:u, :v])
    end
    # Call the standard inner update function for the velocity solve
    invoke(inner_update!, Tuple{AbstractModel}, model)
    # Sync rheology fields used near rank interfaces, then rebuild operators with correct halo data.
    halo_exchange!(model; fields=[:β, :βeff, :ηav])
    update_rheological_operators!(model)
    return model
end

function core_inner_masks(model::Model{<:Any, <:Any, <:MPISpec})
    @unpack gu, gv = model.fields
    th, rh, bh, lh = get_halos(model.spec)

    u_core_mask = falses(size(gu.mask_inner))
    u_core_mask[(1+lh):(size(u_core_mask, 1)-rh), (1+th):(size(u_core_mask, 2)-bh)] .= true
    u_core_inner = u_core_mask[gu.mask_inner]

    v_core_mask = falses(size(gv.mask_inner))
    v_core_mask[(1+lh):(size(v_core_mask, 1)-rh), (1+th):(size(v_core_mask, 2)-bh)] .= true
    v_core_inner = v_core_mask[gv.mask_inner]

    return u_core_inner, v_core_inner
end

"""
    precondition!(model::Model{<:Any, <:Any, <:MPISpec})

Solves the linear system using an iterative overlapping Schwarz method across the distributed domain.

**Iterate** up to `niterations` times:
*   **Local Solve**: Applies the standard `precondition!` locally on each rank.
*   **Interface Update**: Exchanges updated velocities with neighbouring ranks. 
    Depending on the model configuration, this uses either Additive Schwarz with Partition-of-Unity 
    (`pou = true`) or standard Restricted Additive Schwarz (`pou = false`).
*   **Check Convergence**: Computes the global relative residual at the end of each iteration.
    To avoid double-counting, only the non-overlapping "core" unknowns are used in the global MPI reduction.
    *   Each rank computes the squared norm of its local residual contribution using **core unknowns only**
        (excluding overlap/halo regions) to prevent double-counting.
    *   Values are aggregated across all processes using `MPI.Allreduce` with `MPI.SUM`.
    *   Solver exits early if the global relative residual meets the Picard tolerance.
"""
function precondition!(model::Model{<:Any, <:Any, <:MPISpec})
    @unpack niterations, pou = model.spec
    @unpack solver_params = model

    converged = false
    global_rel_resid = Inf

    for iteration = 1:niterations
        if (iteration > 1) && (model.spec.rank == 0)
            @debug "Schwarz iteration $iteration"
        end

        if pou
            # Snapshot only needed when damping mixes in the pre-solve velocity.
            if iszero(model.spec.damping)
                local_precondition!(model)
                mpi_pou_weighted_prolong_velocities!(
                    model,
                    model.fields.gu.u,  # unused when damping == 0
                    model.fields.gv.v,
                )
            else
                u0 = copy(model.fields.gu.u)
                v0 = copy(model.fields.gv.v)
                local_precondition!(model)
                mpi_pou_weighted_prolong_velocities!(model, u0, v0)
            end
        else
            local_precondition!(model)
            halo_exchange!(model; fields=[:u, :v])
        end

        # Check convergence after each iteration
        x = get_start_guess(model)
        op = get_op(model)
        b = get_rhs(model)
        resid = get_resid(x, op, b)

        # Global Residual Check (core-only):
        # exclude overlap halos so each physical unknown is counted once globally.
        @unpack gu, gv = model.fields
        u_core_inner, v_core_inner = core_inner_masks(model)

        # Calculate squared norms locally on core unknowns only
        local_resid_sq = sum(abs2, @view resid[1:gu.ni][u_core_inner]) +
                         sum(abs2, @view resid[(gu.ni+1):(gu.ni+gv.ni)][v_core_inner])
        local_rhs_sq = sum(abs2, @view b[1:gu.ni][u_core_inner]) +
                       sum(abs2, @view b[(gu.ni+1):(gu.ni+gv.ni)][v_core_inner])

        # Sum squared norms across all ranks
        global_resid_sq = MPI.Allreduce(local_resid_sq, MPI.SUM, model.spec.comm)
        global_rhs_sq = MPI.Allreduce(local_rhs_sq, MPI.SUM, model.spec.comm)

        # Compute global relative residual (guard against degenerate zero-RHS case)
        global_rel_resid = iszero(global_rhs_sq) ? sqrt(global_resid_sq) : sqrt(global_resid_sq) / sqrt(global_rhs_sq)

        converged = global_rel_resid < solver_params.tol_picard
        if converged
            if model.spec.rank == 0
                # Early exit avoids remaining PoU/RAS prolongs this Picard step.
                @debug "Schwarz converged early at iteration $iteration / $niterations (rel_resid=$global_rel_resid)"
            end
            break
        end
    end

    if model.spec.rank == 0
        if converged
            @debug "Picard Check: Schwarz early-exit OK; Global Relative Residual = $global_rel_resid (Tol = $(solver_params.tol_picard))"
        else
            @debug "Picard Check: used all $niterations Schwarz iterations; Global Relative Residual = $global_rel_resid (Tol = $(solver_params.tol_picard))"
        end
    end

    return converged, global_rel_resid
end
