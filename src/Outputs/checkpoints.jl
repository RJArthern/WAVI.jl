using JLD2

using WAVI.Parameters: TimesteppingParams
using WAVI.Time: Clock
using WAVI: AbstractModel, AbstractSpec
using WAVI.Fields: InitialConditions, initial_conditions_from_fields, copy_initial_conditions_to_fields!
using WAVI.Processes: update_surface_elevation!, update_geometry_on_uv_grids!,
    update_height_above_floatation!, update_grounded_fraction_on_huv_grids!,
    update_shelf_basal_melt!, update_glen_b!,
    update_dsdh!, update_basal_hydrology!, update_av_viscosity!, update_quadrature_falpha!,
    update_βeff!, update_βeff_on_uv_grids!, update_rheological_operators!,
    update_velocities_on_h_grid!, update_shelf_strain_rate!, update_surf_speed!,
    update_model_wavelets!
using WAVI.SurfaceMassBalance: update_accumulation_rate!

"""
    CheckpointGridGeometry

Size, spacing, and origin of the mesh, stored so pickup can check that the
live model matches the checkpoint. Not a full `Grid`.
"""
struct CheckpointGridGeometry{N <: Integer, T <: Real}
    nx::N
    ny::N
    nσ::N
    dx::T
    dy::T
    x0::T
    y0::T
end

"""Ice-state snapshot (clock, grid, bed, IC arrays, basal drag β). Not a live `Model`."""
struct CheckpointState{C, G, B, I, βT}
    clock::C
    grid::G
    bed::B
    initial_conditions::I
    beta::βT
end

grid_geometry(grid) = CheckpointGridGeometry(grid.nx, grid.ny, grid.nσ, grid.dx, grid.dy, grid.x0, grid.y0)

function dict_get(d, name::Symbol)
    haskey(d, name) && return d[name]
    haskey(d, string(name)) && return d[string(name)]
    return nothing
end

function checkpoint_ics_dict(ics)
    return Dict{String, Any}(string(name) => getfield(ics, name) for name in fieldnames(typeof(ics)))
end

function snapshot_clock(clock::Clock)
    return Dict{String, Any}("n_iter" => clock.n_iter, "time" => clock.time, "ref_time" => clock.ref_time)
end

function snapshot_clock(d::AbstractDict)
    kwargs = Pair{Symbol, Any}[]
    for name in (:n_iter, :time, :ref_time)
        value = dict_get(d, name)
        value === nothing || push!(kwargs, name => value)
    end
    dict_get(d, :n_iter) === nothing && error("Checkpoint clock is missing n_iter")
    dict_get(d, :time) === nothing && error("Checkpoint clock is missing time")
    return Clock(; kwargs...)
end

function snapshot_grid(grid::CheckpointGridGeometry)
    return Dict{String, Any}(
        "nx" => grid.nx, "ny" => grid.ny, "nσ" => grid.nσ,
        "dx" => grid.dx, "dy" => grid.dy, "x0" => grid.x0, "y0" => grid.y0,
    )
end

function snapshot_grid(d::AbstractDict)
    return CheckpointGridGeometry(d["nx"], d["ny"], d["nσ"], d["dx"], d["dy"], d["x0"], d["y0"])
end

snapshot_clock_value(clock::Clock) = clock
snapshot_clock_value(d::AbstractDict) = snapshot_clock(d)
snapshot_grid_value(grid::CheckpointGridGeometry) = grid
snapshot_grid_value(d::AbstractDict) = snapshot_grid(d)

function checkpoint_state(fields, grid, bed, clock::Clock)
    ics = initial_conditions_from_fields(fields)
    return CheckpointState(
        snapshot_clock(clock),
        snapshot_grid(grid_geometry(grid)),
        copy(bed),
        checkpoint_ics_dict(ics),
        copy(fields.gh.β),
    )
end

function copy_snapshot_into_fields!(fields, state::CheckpointState)
    d = state.initial_conditions
    d isa AbstractDict || (d = checkpoint_ics_dict(d))
    present = [name for name in fieldnames(InitialConditions) if dict_get(d, name) !== nothing]
    :initial_thickness in present || error("Checkpoint is missing initial_thickness")
    kwargs = Pair{Symbol, Any}[name => dict_get(d, name) for name in present]
    copy_initial_conditions_to_fields!(fields, InitialConditions(; kwargs...), present)
    if size(state.beta) != size(fields.gh.β)
        error("Checkpoint β has size $(size(state.beta)), live field has $(size(fields.gh.β))")
    end
    fields.gh.β .= state.beta
    return nothing
end

function assert_checkpoint_matches!(grid, bed, state::CheckpointState)
    expected = grid_geometry(grid)
    got = snapshot_grid_value(state.grid)
    if got != expected
        error("Checkpoint grid does not match the live model.")
    end
    if size(state.bed) != size(bed) || state.bed != bed
        error("Checkpoint bed does not match the live model.")
    end
    return nothing
end

function restore_checkpoint_state!(model::AbstractModel, state::CheckpointState)
    assert_checkpoint_matches!(model.grid, model.fields.gh.b, state)
    copy_snapshot_into_fields!(model.fields, state)
    return model
end

"""
Rebuild geometry and derived operators after arrays have been copied in.

"""
function finalise_checkpoint_restore!(model::AbstractModel, clock::Clock)
    update_surface_elevation!(model)
    update_geometry_on_uv_grids!(model)
    update_height_above_floatation!(model)
    update_grounded_fraction_on_huv_grids!(model)
    update_accumulation_rate!(model, clock)
    update_shelf_basal_melt!(model, clock)
    update_glen_b!(model)
    update_dsdh!(model)
    update_basal_hydrology!(model; update_basal_water_thickness = false)
    update_av_viscosity!(model)
    update_quadrature_falpha!(model)
    update_βeff!(model)
    update_βeff_on_uv_grids!(model)
    update_rheological_operators!(model)
    update_velocities_on_h_grid!(model)
    update_shelf_strain_rate!(model)
    update_surf_speed!(model)
    update_model_wavelets!(model)
    return model
end

function checkpoint_path(timestepping_params::TimesteppingParams, output_params::OutputParams)
    if timestepping_params.chkpt_path != "./"
        return timestepping_params.chkpt_path
    elseif output_params.output_path != "./"
        return output_params.output_path
    else
        return timestepping_params.chkpt_path
    end
end

checkpoint_filename(n_iter::Integer) = string("Chkpt_", lpad(n_iter, 10, "0"), ".jld2")

function should_write_checkpoint(timestepping_params::TimesteppingParams, clock::Clock)
    return timestepping_params.chkpt_freq != Inf &&
           clock.n_iter > 0 &&
           mod(clock.n_iter, timestepping_params.n_iter_chkpt) == 0
end

function write_checkpoint!(model::AbstractModel, timestepping_params::TimesteppingParams, output_params::OutputParams, clock::Clock)
    return write_checkpoint!(model.spec, model, timestepping_params, output_params, clock)
end

function write_checkpoint!(
    ::AbstractSpec,
    model,
    timestepping_params::TimesteppingParams,
    output_params::OutputParams,
    clock::Clock,
)
    state = checkpoint_state(model.fields, model.grid, model.fields.gh.b, clock)
    save_collected_checkpoint(state, timestepping_params, output_params, clock)
    return nothing
end

function save_collected_checkpoint(
    state::CheckpointState,
    timestepping_params::TimesteppingParams,
    output_params::OutputParams,
    clock::Clock,
)
    path = checkpoint_path(timestepping_params, output_params)
    isdir(path) || mkpath(path)
    fname = joinpath(path, checkpoint_filename(clock.n_iter))
    @save fname checkpoint_state = state
    @info "Permanent checkpoint at timestep number $(clock.n_iter) - $(fname)"
    return nothing
end

function load_checkpoint_state(
    timestepping_params::TimesteppingParams,
    output_params::OutputParams;
    allow_legacy::Bool = true,
)
    path = checkpoint_path(timestepping_params, output_params)
    fname = joinpath(path, checkpoint_filename(timestepping_params.niter0))
    isfile(fname) || throw_missing_checkpoint(fname, path)
    dict = JLD2.load(fname)
    return checkpoint_state_from_file(dict, fname, timestepping_params.niter0; allow_legacy = allow_legacy)
end

"""
    normalise_checkpoint_state(state, fname, n_iter)

Turn clock and grid into typed values, and ice arrays into a dict, so restore
always sees the same shape of `CheckpointState`. Warns if the file's iteration
does not match `niter0`.
"""
function normalise_checkpoint_state(state::CheckpointState, fname, n_iter)
    clock = snapshot_clock_value(state.clock)
    grid = snapshot_grid_value(state.grid)
    ics = state.initial_conditions
    ics isa AbstractDict || (ics = checkpoint_ics_dict(ics))
    if clock.n_iter != n_iter
        @warn "Checkpoint iteration $(clock.n_iter) does not match niter0=$n_iter ($(fname))"
    end
    return CheckpointState(clock, grid, state.bed, ics, state.beta)
end

function checkpoint_state_from_file(dict, fname, n_iter; allow_legacy::Bool)
    if haskey(dict, "checkpoint_state")
        state = dict["checkpoint_state"]
        state isa CheckpointState || error("Invalid checkpoint_state in $(fname)")
        return normalise_checkpoint_state(state, fname, n_iter)
    end

    if haskey(dict, "model") || haskey(dict, "simulation")
        allow_legacy || error(
            "Checkpoint $(fname) stores a serialised Model. " *
            "This pickup requires a snapshot checkpoint, not a Model-in-JLD2 file.",
        )
        return checkpoint_state_from_legacy(dict, fname, n_iter)
    end

    error("Unrecognised checkpoint format: $(fname)")
end

function checkpoint_state_from_legacy(dict, fname, n_iter)
    if haskey(dict, "simulation")
        model = dict["simulation"].model
        clock = dict["simulation"].clock
    else
        model = dict["model"]
        clock = dict["clock"]
    end
    return normalise_checkpoint_state(
        checkpoint_state(model.fields, model.grid, model.fields.gh.b, clock),
        fname,
        n_iter,
    )
end

function throw_missing_checkpoint(fname, path)
    if isdir(path)
        rank_files = filter(name -> occursin(r"Chkpt_.*_Rank", name), readdir(path))
        if !isempty(rank_files)
            error(
                "Checkpoint file $(fname) not found. " *
                "Per-rank checkpoint files are not supported. Found $(rank_files[1]).",
            )
        end
    end
    error("Checkpoint file $(fname) not found")
end
