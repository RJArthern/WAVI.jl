mutable struct Simulation{T<:Real, N<:Int,TS,O,C} 
    model::Model{T,N}
    timestepping_params::TS
    output_params::O
    clock::C
end

function Simulation(;
                    model = nothing,
                    timestepping_params = nothing,
                    output_params = OutputParams(),
                    pickup_model_update_flag = false,
                    pickup_output_update_flag = false)

    (timestepping_params !== nothing) || throw(ArgumentError("You must specify a timestepping parameters"))

    #compute number of timesteps per output (should be robust for Inf output frequency)
    output_params = set_n_iter_out!(output_params,timestepping_params.dt, timestepping_params.n_iter_total)

    #initialize the clock
    clock = Clock(n_iter = 0, time = 0.0)

    #set the timestep in model parameters (fudge to allow model to see the timestep in velocity solve)
    model = set_dt_in_model!(model, timestepping_params.dt)

    #build the simulation
    simulation = Simulation(model, timestepping_params, output_params, clock)

    #pickup?
    pickup!(simulation, pickup_model_update_flag, pickup_output_update_flag)

    return simulation 

  
    
end
f() =1 
include("run_simulation.jl")

function set_dt_in_model!(model, dt)
    model = @set model.params.dt = dt
    return model
end


function set_n_iter_out!(output_params, dt,n_iter_total)
    output_params.output_freq == Inf ? n_iter_out = (n_iter_total + 1) : n_iter_out = round(Int,output_params.output_freq/dt)
    output_params = @set output_params.n_iter_out = n_iter_out
    return output_params
end

function pickup!(simulation, pickup_model_update_flag, pickup_output_update_flag)
    @unpack timestepping_params = simulation
    if timestepping_params.niter0 > 0
        n_iter_string =  lpad(timestepping_params.niter0, 10, "0"); #filename as a string with 10 digits
        println("detected niter0 > 0 (niter0 = $(timestepping_params.niter0)). Looking for pickup...")
        try 
            @load string("PChkpt_",n_iter_string, ".jld2") simulation
            println("Pickup successful")
            #simulation = @set simulation.timestepping_params = timestepping_params
            simulation.timestepping_params = timestepping_params
            println("Updated the timestepping parameters in the simulation...")
        
            #update model and output if the flags are specified
            if pickup_model_update_flag
                simulation = @set simulation.model = model
                simulation.model = model
                println("WARNING: model updated in simulation after pickup")
            end
            if pickup_output_update_flag
                simulation = @set simulation.output_params = output_params
                simulation.output_params = output_params
                println("WARNING: output parameters updated in simulation after pickup")
            end
        catch 
            Throw(error("Pickup error, terminating run"))
        end
    end
    return simulation
end
    

