# MPISpec: gather snapshot fields to root, save one portable file, scatter on pickup.
# Dispatch is on write_checkpoint!(::MPISpec) and pickup_model(::MPISpec).

using MPI

using WAVI.Outputs: OutputParams, checkpoint_state, save_collected_checkpoint,
    load_checkpoint_state, copy_snapshot_into_fields!, assert_checkpoint_matches!,
    finalise_checkpoint_restore!
import WAVI.Fields: CHECKPOINT_FIELD_PATHS, gridfield_get
using WAVI.Time: Clock
using WAVI.Parameters: TimesteppingParams

checkpoint_global_path(field_path::Tuple) = Symbol[:global_fields, field_path...]

function gather_checkpoint_fields!(model::AbstractModel{<:Any, <:Any, <:MPISpec})
    for path in CHECKPOINT_FIELD_PATHS
        collect_mpi_field!(model, checkpoint_global_path(path))
    end
    return nothing
end

function scatter_checkpoint_fields!(model::AbstractModel{<:Any, <:Any, <:MPISpec})
    for path in CHECKPOINT_FIELD_PATHS
        mpi_fill_local_from_global!(
            model,
            checkpoint_global_path(path),
            gridfield_get(model.fields, path),
        )
    end
    return nothing
end

function WAVI.Outputs.write_checkpoint!(
    ::MPISpec,
    model,
    timestepping_params::TimesteppingParams,
    output_params::OutputParams,
    clock::Clock,
)
    gather_checkpoint_fields!(model)
    ok = zeros(Int32, 1)
    if model.spec.rank == 0
        try
            gf = model.spec.global_fields
            state = checkpoint_state(gf, model.spec.global_grid, gf.gh.b, clock)
            save_collected_checkpoint(state, timestepping_params, output_params, clock)
            ok[1] = 1
        catch e
            @error "MPI root failed to write checkpoint" exception = (e, catch_backtrace())
            ok[1] = 0
        end
    end
    MPI.Bcast!(ok, model.spec.comm)
    ok[1] == 1 || error("Checkpoint write failed on MPI root")
    return nothing
end

function WAVI.Simulations.pickup_model(
    ::MPISpec,
    model::AbstractModel,
    timestepping_params::TimesteppingParams,
    output_params::OutputParams,
)
    @info "detected niter0 > 0 (niter0 = $(timestepping_params.niter0)). Looking for pickup..."
    comm = model.spec.comm
    ok = zeros(Int32, 1)
    n_iter_buf = zeros(Int, 1)
    time_buf = zeros(Float64, 2)

    if model.spec.rank == 0
        try
            state = load_checkpoint_state(timestepping_params, output_params; allow_legacy = false)
            assert_checkpoint_matches!(model.spec.global_grid, model.spec.global_fields.gh.b, state)
            copy_snapshot_into_fields!(model.spec.global_fields, state)
            n_iter_buf[1] = state.clock.n_iter
            time_buf[1] = state.clock.time
            time_buf[2] = state.clock.ref_time
            ok[1] = 1
        catch e
            @error "MPI root failed to load checkpoint" exception = (e, catch_backtrace())
            ok[1] = 0
        end
    end

    MPI.Bcast!(ok, comm)
    ok[1] == 1 || error("Checkpoint pickup failed on MPI root")
    MPI.Bcast!(n_iter_buf, comm)
    MPI.Bcast!(time_buf, comm)
    scatter_checkpoint_fields!(model)
    clock = Clock(n_iter = n_iter_buf[1], time = time_buf[1], ref_time = time_buf[2])
    ok[1] = 1
    try
        climate_clock = WAVI.Simulations.last_completed_climate_forcing_clock(clock, timestepping_params)
        WAVI.Simulations.update_model_climate_forcing!(model, climate_clock)
        finalise_checkpoint_restore!(model, clock)
    catch e
        @error "MPI rank failed to finalise checkpoint pickup" exception = (e, catch_backtrace())
        ok[1] = 0
    end
    ok[1] = MPI.Allreduce(ok[1], MPI.MIN, comm)
    ok[1] == 1 || error("Checkpoint pickup failed during restore")
    model.spec.rank == 0 && println("Pickup successful")
    return (model, clock)
end
