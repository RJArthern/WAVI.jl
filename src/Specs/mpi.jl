export MPISpec

using JLD2
using MPI
using Parameters

using WAVI.Parameters

using WAVI.Architectures: AbstractArchitecture, CPU, GPU

import WAVI: AbstractGrid, AbstractMeltRate, AbstractSurfaceMassBalance, AbstractFracture, AbstractSlidingLaw,
                   AbstractBasalHydrology, AbstractThermoDynamics, AbstractModel
import WAVI.Architectures: architecture, child_architecture, assign_local_device!, on_architecture
import WAVI.Deferred: Collector, register_item!, field_extractor
import WAVI.Fields: GridField, InitialConditions
import WAVI.Grids: Grid, reconstruct_on_grid, reconstruct_on_subdomain
import WAVI.MeltRates: UniformMeltRate
import WAVI.Models: Model, get_bed_elevation
import WAVI.Outputs: write_outputs, zip_output, OutputParams, checkpoint_filename, load_checkpoint, write_checkpoint!,
                    checkpoint_path, with_cleared_stencil_scratch, is_output_step
import WAVI.Parameters: TimesteppingParams
import WAVI.Processes: update_state!, update_model_velocities!, update_velocities!, inner_update!, inner_update_fields!,
                    precondition!, update_preconditioner!, update_rheological_operators!,
                    get_start_guess, get_op, get_rhs
import WAVI.Utilities: stencil_scratch!, get_resid!, _host, _is_host_array
import WAVI.Simulations: run_simulation!, timestep!, update_model_climate_forcing!
import WAVI.Time: Clock

# FIXME: important to realise that this specification has become a complex structure to house many things that should be baked into the model structurally
#  not least the global grid and fields 

"""
Reusable PoU weights, velocity workspaces, and strip packs for neighbour prolong.
Allocated once per local velocity shape; reused across Schwarz iterations.

Weights and work arrays follow the local velocity backend. Send/recv strips are
host `Vector`s on purpose. Device contrib packs into `dev_halo_strip` then
`copyto!` only that strip onto the host buffers for MPI. This is not CUDA-aware MPI.
"""
mutable struct MPIPoUScratch{T <: AbstractFloat}
    ωu::AbstractMatrix{T}
    ωv::AbstractMatrix{T}
    work_u::AbstractMatrix{T}
    work_v::AbstractMatrix{T}
    send_l::Vector{T}
    recv_l::Vector{T}
    send_r::Vector{T}
    recv_r::Vector{T}
    send_t::Vector{T}
    recv_t::Vector{T}
    send_b::Vector{T}
    recv_b::Vector{T}
    # GPU pack/unpack workspace for one PoU overlap edge (`nothing` until a GPU prolong).
    dev_halo_strip::Any
end

function MPIPoUScratch(ωu::AbstractMatrix{T}, ωv::AbstractMatrix{T}) where {T <: AbstractFloat}
    empty = T[]
    return MPIPoUScratch{T}(
        ωu, ωv, similar(ωu), similar(ωv),
        copy(empty), copy(empty), copy(empty), copy(empty),
        copy(empty), copy(empty), copy(empty), copy(empty),
        nothing,
    )
end

"""
Reusable RAS halo packs, blend weights, and L0 snapshots.
Allocated on first `halo_exchange!`; send/recv strips grow to the largest field.

Packs, L0 snapshots, and blend weights are host `Array`s on purpose.
Device fields pack into `dev_halo_strip` then `copyto!` the halo strip onto these host
buffers for MPI. This is not CUDA-aware MPI.
"""
mutable struct MPIHaloScratch{T <: AbstractFloat}
    send_l::Vector{T}
    recv_l::Vector{T}
    send_r::Vector{T}
    recv_r::Vector{T}
    send_t::Vector{T}
    recv_t::Vector{T}
    send_b::Vector{T}
    recv_b::Vector{T}
    l0_h::Matrix{T}
    l0_u::Matrix{T}
    l0_v::Matrix{T}
    W_left::Vector{T}
    W_right::Vector{T}
    W_top::Matrix{T}
    W_bottom::Matrix{T}
    # GPU pack/unpack workspace for one halo edge (`nothing` until a GPU field is exchanged).
    dev_halo_strip::Any
end

function MPIHaloScratch(::Type{T}, halo::Integer, damping) where {T <: AbstractFloat}
    empty = T[]
    d = T(damping)
    return MPIHaloScratch{T}(
        copy(empty), copy(empty), copy(empty), copy(empty),
        copy(empty), copy(empty), copy(empty), copy(empty),
        zeros(T, 0, 0), zeros(T, 0, 0), zeros(T, 0, 0),
        fill(d, halo), fill(d, halo), fill(d, 1, halo), fill(d, 1, halo),
        nothing,
    )
end

"""
Struct to represent the MPI parallel specification of a model.

Fields:
    px, py: Number of processes in x and y directions
    halo: Number of halo cells on each side of the subgrid
    rank: MPI rank of the current process
    comm: MPI communicator for the current process
    coords: Cartesian coordinates of the current process as (x, y)
    top, right, bottom, left: Neighbouring processes in the x and y directions; MPI_PROC_NULL if none
    global_size: Total number of processes
    global_comm: MPI communicator for the global grid
    global_grid: Global grid to be used for the model
    damping: Damping factor for halo exchange (default=0.0)
    niterations: Number of Schwarz iterations per Picard iteration (default=5)
    field_collector: Field collector for the model
    pou_scratch: Cached PoU weights (device copy on GPU ranks) and host MPI strip buffers (filled on first prolong)
    halo_scratch: Cached RAS host halo packs plus a device halo strip workspace (filled on first halo_exchange!)
    core_inner: Cached core-only inner masks for the global residual (host on CPU, device copy on GPU)
    child_architecture: Where each rank's local arrays live (`CPU()` or `GPU()`)
    local_spec: Optional intraprocess ThreadedSpec for the rank-local solve (default=nothing to wavelet)
"""
mutable struct MPISpec{N <: Integer, T <: Number, M <: MPI.Comm, G <: AbstractGrid, A <: AbstractArchitecture} <: AbstractDecompSpec
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
    coords::NTuple{2,N}

    # Neighbourhood information (MPI_PROC_NULL if no neighbour)
    top::N
    right::N
    bottom::N
    left::N
    pou::Bool
    damping::T
    niterations::N  # Number of Schwarz iterations per Picard iteration
    pou_scratch::Union{Nothing, MPIPoUScratch{T}}
    halo_scratch::Union{Nothing, MPIHaloScratch{T}}
    core_inner::Union{Nothing, Tuple{AbstractVector{Bool}, AbstractVector{Bool}}}
    child_architecture::A
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
    child_architecture: `CPU()` (default) or `GPU()` after `using CUDA`. One GPU per rank.
        Do not combine `GPU()` with `local_spec = ThreadedSpec(...)`.
    """ function MPISpec(
        px::Integer,
        py::Integer,
        halo::Integer,
        grid::G;
        pou::Bool = true,
        damping::AbstractFloat = 0.0,
        niterations::Integer = 5,
        child_architecture::AbstractArchitecture = CPU(),
        local_spec::Union{Nothing, ThreadedSpec} = nothing,
    ) where {G <: AbstractGrid}
        (px < 1 || py < 1 || halo < 0) &&
            throw(ArgumentError("Invalid parameters specified for MPISpec"))

        MPI.Initialized() || MPI.Init()
        comm = MPI.COMM_WORLD
        rank = MPI.Comm_rank(comm)
        size = MPI.Comm_size(comm)
        @debug "Creating dimensions of $(size) with ($(px), $(py))"

        if child_architecture isa GPU && local_spec !== nothing
            throw(ArgumentError(
                "ThreadedSpec cannot be used as local_spec when child_architecture is a GPU. " *
                "Use local_spec = nothing (the default) for GPU ranks.",
            ))
        end

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

        if !(child_architecture isa CPU)
            node_comm = MPI.Comm_split_type(comm, MPI.COMM_TYPE_SHARED, rank)
            assign_local_device!(
                child_architecture,
                MPI.Comm_rank(node_comm),
                MPI.Comm_size(node_comm),
            )
        end

        return new{Int, Float64, MPI.Comm, G, typeof(child_architecture)}(
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
            (Int(x_coord), Int(y_coord)),
            Int(top),
            Int(right),
            Int(bottom),
            Int(left),
            pou,
            damping,
            niterations,
            nothing,  # pou_scratch filled on first PoU prolong
            nothing,  # halo_scratch filled on first RAS exchange
            nothing,  # core_inner filled on first Schwarz residual
            child_architecture,
            local_spec,
        )
    end
end

architecture(spec::MPISpec) = spec.child_architecture
child_architecture(spec::MPISpec) = spec.child_architecture

include("MPI/utils.jl")
include("MPI/exchanges.jl")
include("MPI/outputs.jl")
include("MPI/mpi_checkpoints.jl")

"""
    mpi_allocate_global_fields!(spec, grid, bed_array; kwargs...)

Create the full-domain arrays that rank 0 uses when writing output.

Only rank 0 stores `spec.global_fields`. Other ranks leave it as `nothing` and send
their patch when fields are gathered. These arrays are only for output, they are not
used for the velocity solve or other physics.
"""
function mpi_allocate_global_fields!(
    spec::MPISpec,
    grid::AbstractGrid,
    bed_array;
    initial_conditions::InitialConditions = InitialConditions(),
    params::Params = Params(),
    solver_params::SolverParams = SolverParams(),
)
    spec.rank == 0 || return nothing
    # Host arrays only: collect/output stay on rank 0 CPU even when local fields are on a GPU.
    spec.global_fields = GridField(
        grid,
        bed_array;
        initial_conditions,
        params,
        solver_params,
        mpi_rank = zeros(Float64, grid.nx, grid.ny),
        assembly_buffer = true,
    )
    return spec.global_fields
end

"""
    mpi_restore_global_fields!(spec, saved_global_fields, model)

Put rank 0's output arrays back after loading a checkpoint.

Loading a checkpoint builds a new `MPISpec`, so those arrays would otherwise be
missing. Rank 0 reuses the saved arrays when they are in the file, or allocates
empty ones. Other ranks do nothing.
"""
function mpi_restore_global_fields!(spec::MPISpec, saved_global_fields, model)
    spec.rank == 0 || return nothing
    if saved_global_fields isa GridField
        spec.global_fields = saved_global_fields
        return spec.global_fields
    end
    grid = spec.global_grid
    return mpi_allocate_global_fields!(
        spec,
        grid,
        zeros(grid.nx, grid.ny);
        params = model.params,
        solver_params = model.solver_params,
    )
end

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

    bounds = (x_start, x_end, y_start, y_end)
    
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

    # Slice global arrays first, then expand remaining scalars onto the local grid.
    #expand scalar paramaters onto grid
    # dt cannot be copied via the external constructor so we create the structure directly
    local_params = reconstruct_on_grid(reconstruct_on_subdomain(params, grid, bounds), local_grid)
    #Replace all NaN entries with defaults from params on correct grid
    #trim initial conditions to local domain
    local_initial_conditions = reconstruct_on_grid(reconstruct_on_subdomain(initial_conditions, grid, bounds), local_params, local_grid)
    #expand spatial parameters onto grid
    local_shelf_melt_rate = reconstruct_on_grid(reconstruct_on_subdomain(shelf_melt_rate, grid, bounds), local_grid)
    local_surface_mass_balance = reconstruct_on_grid(reconstruct_on_subdomain(surface_mass_balance, grid, bounds), local_grid)
    local_fracture = reconstruct_on_grid(reconstruct_on_subdomain(fracture, grid, bounds), local_grid)
    local_sliding_law = reconstruct_on_grid(reconstruct_on_subdomain(sliding_law, grid, bounds), local_grid)
    local_basal_hydrology = reconstruct_on_grid(reconstruct_on_subdomain(basal_hydrology, grid, bounds), local_grid)
    local_thermo_dynamics = reconstruct_on_grid(reconstruct_on_subdomain(thermo_dynamics, grid, bounds), local_grid)

    if typeof(bed_elevation) <: AbstractArray
        bed_array = bed_elevation[x_start:x_end, y_start:y_end]
    else
        bed_array = get_bed_elevation(bed_elevation, local_grid)
    end

    # Create local mpi_rank field filled with this rank's number
    local_mpi_rank = fill(Float64(rank), nx_local, ny_local)

    fields = GridField(local_grid, bed_array; initial_conditions=local_initial_conditions, params=local_params, solver_params, mpi_rank=local_mpi_rank)
    arch = architecture(spec)
    if !(arch isa CPU)
        # Rank-0 global_fields stay on the host for collect/output (see mpi_allocate_global_fields!).
        fields = on_architecture(arch, fields)
    end
    model = Model(local_grid, fields, local_params, solver_params, spec, local_shelf_melt_rate, local_surface_mass_balance, local_fracture, local_sliding_law, local_basal_hydrology, 
    local_thermo_dynamics, verbose)

    # Assembly buffers stay on rank 0 only. Other ranks send patches during collect.
    # Global assembly buffers for gather/Bcast only: sized arrays, no IC/operator cost.
    # (Physics lives on each rank's local GridField; bed is retained for static outputs.)
    if rank == 0
        global_bed = typeof(bed_elevation) <: AbstractArray ? bed_elevation : get_bed_elevation(bed_elevation, grid)
        mpi_allocate_global_fields!(spec, grid, global_bed; initial_conditions, params, solver_params)
    end

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

    # Collect AFTER thickness update to match BasicSpec output timing.
    # Gather registered global fields only when write_output will actually write.
    if is_output_step(output_params, clock)
        collect!(model.spec.field_collector, model)
    end
    write_outputs(model, timestepping_params, output_params, clock)
    if is_output_step(output_params, clock)
        clear!(model.spec.field_collector)
    end
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
        return collect_mpi_field!(model, path)
    end
    extractor = field_extractor(join(string.(path), "."), accessor, path)
    register_item!(collector, extractor)
end

"""
    inner_update!(model::Model{<:Any, <:Any, <:MPISpec})

Overload for the inner update of the velocity solve.

One-rank (no neighbours) uses the same local path as GPUSpec. With neighbours,
sync halos so viscosity and other calculations have boundary data from
neighbouring ranks. `update_rheological_operators!` runs after that halo
sync because β, βeff, and ηav must include cross-rank values.
"""
function inner_update!(model::Model{<:Any, <:Any, <:MPISpec})
    mpi_has_neighbours(model.spec) || return _inner_update_local!(model)
    return _inner_update_mpi_neighbours!(model)
end

function _inner_update_local!(model::Model{<:Any, <:Any, <:MPISpec})
    inner_update_fields!(model)
    update_rheological_operators!(model)
    return model
end

function _inner_update_mpi_neighbours!(model::Model{<:Any, <:Any, <:MPISpec})
    # RAS ghost sync for velocities; AS-PoU keeps overlap values from the last prolongation.
    if !model.spec.pou
        halo_exchange!(model; fields=[:u, :v])
    end
    inner_update_fields!(model)
    halo_exchange!(model; fields=[:β, :βeff, :ηav])
    update_rheological_operators!(model)
    return model
end

function core_inner_masks(model::Model{<:Any, <:Any, <:MPISpec})
    cached = model.spec.core_inner
    cached !== nothing && return cached

    @unpack gu, gv = model.fields
    th, rh, bh, lh = get_halos(model.spec)

    # Host copies: `mask_inner` lives on the child architecture (GPU ranks).
    mask_u = _host(gu.mask_inner)
    mask_v = _host(gv.mask_inner)

    u_core_mask = falses(size(mask_u))
    u_core_mask[(1+lh):(size(u_core_mask, 1)-rh), (1+th):(size(u_core_mask, 2)-bh)] .= true
    u_core_inner = Vector{Bool}(u_core_mask[mask_u])

    v_core_mask = falses(size(mask_v))
    v_core_mask[(1+lh):(size(v_core_mask, 1)-rh), (1+th):(size(v_core_mask, 2)-bh)] .= true
    v_core_inner = Vector{Bool}(v_core_mask[mask_v])

    proto = gu.u
    if _is_host_array(proto)
        model.spec.core_inner = (u_core_inner, v_core_inner)
    else
        u_dev = similar(proto, Bool, length(u_core_inner))
        v_dev = similar(proto, Bool, length(v_core_inner))
        copyto!(u_dev, u_core_inner)
        copyto!(v_dev, v_core_inner)
        model.spec.core_inner = (u_dev, v_dev)
    end
    return model.spec.core_inner
end

"""
    masked_sum_abs2(x, r, mask)

Sum of squares of `x[r]` where `mask` is true. Host arrays use boolean views.
Device arrays reduce on the same backend; only the scalar returns to the host.
"""
function masked_sum_abs2(x, r, mask)
    if _is_host_array(x)
        return sum(abs2, @view x[r][mask])
    end
    # Dense slice: a view of a device array can scalar-index under boolean mask.
    return sum(abs2, x[r][mask])
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
    model.spec.global_size == 1 && return local_precondition!(model)
    return _schwarz_precondition!(model)
end

function _schwarz_precondition!(model::Model{<:Any, <:Any, <:MPISpec})
    @unpack niterations, pou = model.spec
    @unpack solver_params = model

    converged = false
    global_rel_resid = Inf
    s = stencil_scratch!(model)

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
        resid = s.picard_resid
        get_resid!(resid, x, op, b)

        # Global Residual Check (core-only):
        # exclude overlap halos so each physical unknown is counted once globally.
        # Packed residual may be on GPU. Core masks live on the same backend;
        # only two scalars are Allreduced.
        @unpack gu, gv = model.fields
        u_core_inner, v_core_inner = core_inner_masks(model)
        u_range = 1:gu.ni
        v_range = (gu.ni + 1):(gu.ni + gv.ni)

        # Calculate squared norms locally on core unknowns only
        local_resid_sq = masked_sum_abs2(resid, u_range, u_core_inner) +
                         masked_sum_abs2(resid, v_range, v_core_inner)
        local_rhs_sq = masked_sum_abs2(b, u_range, u_core_inner) +
                       masked_sum_abs2(b, v_range, v_core_inner)

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
