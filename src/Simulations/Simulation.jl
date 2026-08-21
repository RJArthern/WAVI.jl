module Simulations

export Simulation

using Parameters
using Setfield
using ImageFiltering: centered, imfilter, reflect, Fill

using WAVI: AbstractModel, AbstractSimulation, AbstractSpec
using WAVI.Outputs: OutputParams, load_checkpoint_state, restore_checkpoint_state!,
    finalise_checkpoint_restore!
using WAVI.Parameters: TimesteppingParams
using WAVI.Time

struct Simulation{M,TS,O,C} <: AbstractSimulation
    model::M
    timestepping_params::TS
    output_params::O
    clock::C
end

"""
    Simulation(;
            model,
            timestepping_params,
            output_params = OutputParams())

Construct a WAVI.jl Simulation object.

Keyword arguments
=================
- `model`: (required) an instance of a `Model` object. When `timestepping_params.niter0 > 0`,
  restart arrays are copied into this model from the checkpoint; physics objects stay as built.
- `timestepping_params`: (required) an instance of a `TimesteppingParams` object, which stores information relating to timestepping
- `output_params`: an instance of an `OutputParams` object, which stores information relating to outputting of solutions.
  Also used with `timestepping_params` to locate checkpoint files for pickup (see `checkpoint_path` in `Outputs`).
"""
function Simulation(; 
                    model::AbstractModel,
                    timestepping_params::TimesteppingParams,
                    output_params::OutputParams = OutputParams())
    isnothing(timestepping_params) && throw(ArgumentError("You must specify a timestepping parameters"))

    #compute number of timesteps per output (should be robust for Inf output frequency)
    output_params = set_n_iter_out!(output_params, timestepping_params.dt, timestepping_params.n_iter_total)
    #set the timestep in model parameters
    model = set_dt_in_model!(model, timestepping_params.dt)

    if timestepping_params.niter0 > 0
        model, clock = pickup_model(model, timestepping_params, output_params)
    else
        # TODO: is the change from the default relevant - time is now type-variant (Real not Int)
        clock = Clock(n_iter = 0, time = 0.0, ref_time = timestepping_params.ref_time)
    end

    return Simulation(model, timestepping_params, output_params, clock)
end

Simulation(m::AbstractModel, tp::TimesteppingParams; kwargs...) = Simulation(; model=m, timestepping_params=tp, kwargs...)

include("run_simulation.jl")

function set_dt_in_model!(model, dt)
    # TODO: code smell, this should be in Model construction, requires model recreation via SetField and ConstructionBase
    model = @set model.params.dt = dt
    return model
end


function set_n_iter_out!(output_params, dt, n_iter_total)
    # TODO: code smell, this should be in the constructor for OutputParams
    output_params.output_freq == Inf ? n_iter_out = (n_iter_total + 1) : n_iter_out = round(Int, output_params.output_freq/dt)
    output_params = @set output_params.n_iter_out = n_iter_out
    return output_params
end

function pickup_model(model::AbstractModel, timestepping_params::TimesteppingParams, output_params::OutputParams)
    return pickup_model(model.spec, model, timestepping_params, output_params)
end

function pickup_model(
    ::AbstractSpec,
    model::AbstractModel,
    timestepping_params::TimesteppingParams,
    output_params::OutputParams,
)
    @info "detected niter0 > 0 (niter0 = $(timestepping_params.niter0)). Looking for pickup..."
    try
        state = load_checkpoint_state(timestepping_params, output_params)
        restore_checkpoint_state!(model, state)
        update_model_climate_forcing!(model, last_completed_climate_forcing_clock(state.clock, timestepping_params))
        finalise_checkpoint_restore!(model, state.clock)
        println("Pickup successful")
        return (model, state.clock)
    catch e
        @error "Pickup error" exception = (e, catch_backtrace())
        rethrow()
    end
end
    
end
