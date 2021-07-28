struct TimesteppingParams{T <: Real, N <: Integer, TO, C, P}
            niter0 :: N      #starting iteration number
                dt :: T      #timestep
          end_time :: T      #end time of this simulation
                t0 :: T      #start time of this simulation 
        chkpt_freq :: T      #temporary checkpoint frequency
       pchkpt_freq :: T      #permanent checkpoint frequency  
       chkpt_path  :: String #path to location of permanent and tempoary checkpoint output
      n_iter_total :: TO     #total number of timesteps counting from zero
      n_iter_chkpt :: C      #number of iterations per temporary checkpoint
     n_iter_pchkpt :: P      #number of iterations per permanent checkpoint
    step_thickness :: Bool   #toggle whether to step the thickness at each timestep or not (coupling control)
end



"""
TimesteppingParams(;<kwargs>)

Construct a WAVI.jl TimesteppingParams object, which stores parameters relating to timestepping.

Keyword arguments
=================
- `niter0`: Iteration number of the first timestep. niter0 = 0 corresponds to a new simulation, while niter0 > 0 (positive integer) corresponds to a pickup.
- `dt`: simulation timestep
- `end_time`: Simulation termination time
- `t0`: Starting time of the simulation
- `chkpt_freq`: model time interval at which temporary checkpoints are outputted
- `pchkpt_freq`: model time interval at which permanent checkpoints are outputted
- `chkpt_path` : Path to location checkpoint output
- `step_thickness`: Toggle whether to update the ice thickness (true) or not (false) at each timestep
"""
function TimesteppingParams(;
                        niter0 = 0,
                        dt = 1.0,
                        end_time = nothing,
                        n_iter_total = nothing, 
                        t0 = nothing,
                        chkpt_freq = Inf,
                        pchkpt_freq = Inf,
                        chkpt_path = "./",
                        step_thickness = true)


    #initialize t0 (really you should read start time from pickup file)
    t0 = niter0 > 0 ? niter0 * dt : 0 
    t0 = map(typeof(dt), t0)

    #check compatibility of n_iter_total and end_time, and compute them 
    end_time, n_iter_total = compute_iterations_and_end_time(end_time, n_iter_total,dt)

    #compute number of timesteps checkpoint number of timesteps
    chkpt_freq == Inf ? n_iter_chkpt = Inf : n_iter_chkpt  = round(Int, chkpt_freq/dt)
    pchkpt_freq == Inf ? n_iter_pchkpt = Inf : n_iter_pchkpt = round(Int, pchkpt_freq/dt)
    
    #check the output path ends in '/' and exists
    endswith(chkpt_path, "/") || (chkpt_path = string(chkpt_path, "/"))
    if ~isdir(chkpt_path)
        @warn string("Did not find output path ", chkpt_path, ". Any outputs will go to the working directory", pwd())
        chkpt_path = "./"
    end

    return TimesteppingParams(niter0, dt, end_time, t0, chkpt_freq, pchkpt_freq, 
                            chkpt_path,n_iter_total, n_iter_chkpt, n_iter_pchkpt, step_thickness)
end

"""
    compute_iterations_and_end_time(end_time, n_iter_total, dt)

Check input end time and total iterations, and compute those not passed.
"""
function compute_iterations_and_end_time(end_time, n_iter_total, dt)
    (isa(dt, Real) && (dt > 0))  || throw(ArgumentError("timestep dt must be a positive number"))

    if (~(n_iter_total === nothing) && ~(end_time === nothing)) #if both passed, throw error if incompatible
        (end_time ≈ (n_iter_total * dt)) ||  throw(ArgumentError("You have specified both end time (end_time) and total iterations (n_iter_total), but their values are incompatible: end time mustequal n_iter_total * dt"))
    elseif ((n_iter_total === nothing) && ~(end_time === nothing)) #if only end time passed, n_iter_total is the nearest integer
        end_time == Inf ? n_iter_total = Inf : n_iter_total  = round(Int, end_time/dt)
    elseif (~(n_iter_total === nothing) && (end_time === nothing)) #if only number of iterations
        n_iter_total == Inf ? end_time = Inf : end_time  = dt*n_iter_total
    elseif (n_iter_total === nothing) && (end_time === nothing) #if neither is passed
        throw(ArgumentError("You must pass at least one of end_time or n_iter_total"))
    end 

    #check both are positive numbers, now that both end_time and n_iter_total assigned
    (isa(n_iter_total, Integer) && (n_iter_total > 0)) || throw(ArgumentError("n_iter_total must be a positive integer"))
    (isa(end_time, Real) && (end_time > 0))|| throw(ArgumentError("end_time must be a positive number"))
    return end_time, n_iter_total
end