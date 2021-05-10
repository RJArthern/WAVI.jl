struct Simulation{T <: Real,N <: Int,R <: Real} <: AbstractSimulation{T,N,R}
    model::Model{T,N}
    timestepping_params::TimesteppingParams{T,N}
    output_params::OutputParams{T,R}
    clock::Clock{T,N}
end

function Simulation(;
                    model = nothing,
                    timestepping_params = nothing,
                    output_params = Output_params())

    (model !== nothing) || throw(ArgumentError("You must specify a model input"))
    (timestepping_params !== nothing) || throw(ArgumentError("You must specify a model input"))

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
end