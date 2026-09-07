using JLD2
using Logging

using WAVI.Parameters: TimesteppingParams
using WAVI.Time: Clock
using WAVI: AbstractSpec

# Drop JLD2's "stores functions by name" warning; keep the caller's logger as-is.
struct SkipJLD2FunctionWarnings <: AbstractLogger
    parent::AbstractLogger
end
Logging.min_enabled_level(l::SkipJLD2FunctionWarnings) = Logging.min_enabled_level(l.parent)
Logging.catch_exceptions(l::SkipJLD2FunctionWarnings) = Logging.catch_exceptions(l.parent)
Logging.shouldlog(l::SkipJLD2FunctionWarnings, args...) = Logging.shouldlog(l.parent, args...)
function Logging.handle_message(l::SkipJLD2FunctionWarnings, level, message, args...; kwargs...)
    level == Warn && occursin("stores functions by name", string(message)) && return
    Logging.handle_message(l.parent, level, message, args...; kwargs...)
end

"""
    checkpoint_path(timestepping_params, output_params)

Directory for permanent checkpoint files (`Chkpt_*.jld2`).

Uses `timestepping_params.chkpt_path` when set explicitly; otherwise falls back to
`output_params.output_path` so runs that only set an output folder still write and
read checkpoints in the same place.
"""
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

"""
    with_cleared_stencil_scratch(f, model)

Run `f()` with checkpoint-hostile caches omitted from the model.

Temporarily sets `model.fields.stencil_scratch` to `nothing`, and replaces
`model.spec.field_collector` with an empty collector when that field exists.
Scratch holds anonymous matvec / restrict / prolong closures; the MPI field
collector holds anonymous accessors. JLD2 cannot serialise either usefully.
Restore both afterwards so an ongoing run keeps its workspace and registered
outputs.
"""
function with_cleared_stencil_scratch(f, model)
    scratch_ref = model.fields.stencil_scratch
    scratch = scratch_ref[]
    scratch_ref[] = nothing

    spec = model.spec
    collector = nothing
    if hasproperty(spec, :field_collector)
        collector = spec.field_collector
        spec.field_collector = Collector()
    end
    try
        return with_logger(SkipJLD2FunctionWarnings(current_logger())) do
            f()
        end
    finally
        scratch_ref[] = scratch
        if collector !== nothing
            spec.field_collector = collector
        end
    end
end

function write_checkpoint!(model, timestepping_params::TimesteppingParams, output_params::OutputParams, clock::Clock)
    path = checkpoint_path(timestepping_params, output_params)
    if !isdir(path)
        mkpath(path)
    end
    fname = joinpath(path, checkpoint_filename(clock.n_iter))
    with_cleared_stencil_scratch(model) do
        @save fname model=model timestepping_params=timestepping_params clock=clock
    end
    @info "Permanent checkpoint at timestep number $(clock.n_iter) — $(fname)"
    return nothing
end

function load_checkpoint(timestepping_params::TimesteppingParams, output_params::OutputParams)
    path = checkpoint_path(timestepping_params, output_params)
    n_iter = timestepping_params.niter0
    fname = joinpath(path, checkpoint_filename(n_iter))
    dict = load(fname)
    if haskey(dict, "simulation")
        # legacy format (pre #105): @save fname simulation saved everything
        # in one bundle, need to unpack it.
        sim_load = dict["simulation"]
        model = sim_load.model
        clock = sim_load.clock
    else
        # new format (post #105): @save fname model=model clock=clock saved separately
        model = dict["model"]
        clock = dict["clock"]
    end
    if clock.n_iter != n_iter
        @warn "Checkpoint iteration $(clock.n_iter) does not match niter0=$n_iter ($(fname))"
    end
    return model, clock
end


# Default behaviour. So that checkpointing can be overloaded by MPISpec.

checkpoint_filename(spec::AbstractSpec,n_iter::Integer) = checkpoint_filename(n_iter::Integer)

write_checkpoint!(spec::AbstractSpec, model, timestepping_params::TimesteppingParams, output_params::OutputParams, clock::Clock) =
    write_checkpoint!(model, timestepping_params::TimesteppingParams, output_params::OutputParams, clock::Clock)

load_checkpoint(spec::AbstractSpec,timestepping_params::TimesteppingParams, output_params::OutputParams) =
 load_checkpoint(timestepping_params::TimesteppingParams, output_params::OutputParams)

