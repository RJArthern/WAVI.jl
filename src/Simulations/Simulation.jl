struct Simulation{T <: Real,N <: Int,R <: Real} <: AbstractSimulation{T,N,R}
    model::Model{T,N}
    timestepping_params::TimesteppingParams{T,N}
    output_params::OutputParams{T,R}
    clock::Clock{T,N}
end

function Simulation(;
                    model = nothing,
                    timestepping_params = nothing,
                    output_params = OutputParams(),
                    pickup_model_update_flag = false,
                    pickup_output_update_flag = false)

    (timestepping_params !== nothing) || throw(ArgumentError("You must specify a timestepping parameters"))

    if timestepping_params.niter0 > 0
        n_iter_string =  lpad(timestepping_params.niter0, 10, "0"); #filename as a string with 10 digits
        println("detected niter0 > 0 (niter0 = $(timestepping_params.niter0)). Looking for pickup...")
        try 
            @load string("PChkpt_",n_iter_string, ".jld2") simulation
            println("Pickup successful")
            simulation = @set simulation.timestepping_params = timestepping_params
            println("Updated the timestepping parameters in the simulation...")
            
            #update model and output if the flags are specified
            if pickup_model_update_flag
                simulation = @set simulation.model = model
                println("WARNING: model updated in simulation after pickup")
            end
            if pickup_output_update_flag
                simulation = @set simulation.output_params = output_params
                println("WARNING: output parameters updated in simulation after pickup")
            end

            return simulation 
        catch 
            Throw(error("Pickup error, terminating run"))
        end
    
    
    elseif iszero(timestepping_params.niter0)
        println("detected niter0 = 0, initializng clean simulation")
        #initialize the clock
        clock = Clock(n_iter = 0, time = 0.0)

        #set the timestep in model parameters (hack to allow model to see the timestep in velocity solve)
        model = @set model.params.dt = timestepping_params.dt

        #initialize number of steps in output
        if output_params.output_freq !== Inf #set the output number of timesteps, if it has been specifies
            n_iter_out = round(Int, output_params.output_freq/timestepping_params.dt) #compute the output timestep
            output_params = @set output_params.n_iter_out = n_iter_out
        end
        return Simulation(model, timestepping_params, output_params, clock)
    else
        throw(ArgumentError("niter0 must be a non-negative integer"))
    end

    
end

include("run_simulation.jl")