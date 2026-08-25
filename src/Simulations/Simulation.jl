module Simulations

export Simulation

using JLD2
using Parameters
using Setfield
using ImageFiltering: centered, imfilter, reflect, Fill

using WAVI: AbstractModel, AbstractSimulation, AbstractSpec
using WAVI.Models: restore_pickup_architecture!
using WAVI.Outputs: OutputParams, load_checkpoint
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
- `model`: (required) an instance of a `Model` object. When `timestepping_params.niter0 > 0`, the model and clock are
  loaded from the checkpoint file; the passed-in `model` is replaced after pickup.
- `timestepping_params`: (required) an instance of a `TimesteppingParams` object, which stores information relating to timestepping
- `output_params`: an instance of an `OutputParams` object, which stores information relating to outputting of solutions.
  Also used with `timestepping_params` to locate checkpoint files for pickup (see `checkpoint_path` in `Outputs`).
"""
function Simulation(; 
                    model::AbstractModel,
                    timestepping_params::TimesteppingParams,
                    output_params::OutputParams = OutputParams())
    isnothing(timestepping_params) && throw(ArgumentError("You must specify a timestepping parameters"))

    @unpack spec = model

    #compute number of timesteps per output (should be robust for Inf output frequency)
    output_params = set_n_iter_out!(output_params, timestepping_params.dt, timestepping_params.n_iter_total)
    pickup_model, pickup_clock = pickup!(spec, timestepping_params, output_params)

    if ~isnothing(pickup_model)
        spec == pickup_model.spec || error("Model Spec does not match spec from pickup file")
        model, clock = pickup_model, pickup_clock
    else
        # TODO: is the change from the default relevant - time is now type-variant (Real not Int)
        clock = Clock(n_iter = 0, time = 0.0, ref_time = timestepping_params.ref_time)
        #set the timestep in model parameters (fudge to allow model to see the timestep in velocity solve)
        model = set_dt_in_model!(model, timestepping_params.dt)
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

function pickup!(spec::AbstractSpec,timestepping_params::TimesteppingParams, output_params::OutputParams)::Union{Tuple{Model, Clock}, Tuple{Nothing, Nothing}}
    model, clock = nothing, nothing
    if timestepping_params.niter0 > 0
        @info "detected niter0 > 0 (niter0 = $(timestepping_params.niter0)). Looking for pickup..."
        try
            model, clock = load_checkpoint(spec,timestepping_params, output_params)
            model = restore_pickup_architecture!(model)
            println("Pickup successful")
        catch e
            @error "Pickup error: $e"
            error("Pickup error, terminating run")
        end
    end
    return (model, clock)
end
    
end
