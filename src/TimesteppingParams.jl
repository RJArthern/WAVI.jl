
struct TimesteppingParams{T <: Real, N <: Integer, TO, C, P}
                        niter0 :: N      #starting iteration number
                            dt :: T      #timestep
    ntimesteps_velocity_update :: N      #number of substeps at which to update the velocity (i.e. the velocity is updated every dt*ntimesteps_velocity_update)
                      end_time :: T      #end time of this simulation
                            t0 :: T      #start time of this simulation 
                    chkpt_freq :: T      #temporary checkpoint frequency
                   pchkpt_freq :: T      #permanent checkpoint frequency  
                   chkpt_path  :: String #path to location of permanent and tempoary checkpoint output
                  n_iter_total :: TO     #total number of timesteps counting from zero
                  n_iter_chkpt :: C      #number of iterations per temporary checkpoint
                 n_iter_pchkpt :: P      #number of iterations per permanent checkpoint
                step_thickness :: Bool   #toggle whether to step the thickness at each timestep or not (coupling control)
                       verbose :: Bool   #toggle whether or not to output when the timestepping have been performed
end



"""
TimesteppingParams(;
                    niter0 = 0,
                    dt = 1.0,
                    ntimesteps_velocity_update = 1,
                    end_time = 1.0,
                    t0 = nothing,
                    chkpt_freq = Inf,
                    pchkpt_freq = Inf,
                    chkpt_path = './',
                    step_thickness = true,
                    verbose = false)

Construct a WAVI.jl TimesteppingParams object.
TimesteppingParams stores information relating to timestepping.

Keyword arguments
=================
- 'niter0': Iteration number of the first timestep. niter0 = 0 corresponds to a new simulation, while niter0 > 0 (positive integer) corresponds to a pickup.
- 'dt': Model timestep
- 'ntimesteps_velocity_update': number of substeps at which to update the velocity (i.e. the velocity is updated every dt*ntimesteps_velocity_update)
- 'end_time': Simulation termination time
- 't0': Starting time of the simulation
- 'chkpt_freq': Frequency of outputting temporary checkpoints
- 'pchkpt_freq': Frequecy with which permanent checkpoints are pass
- 'chkpt_path' : Path to location checkpoint output
- 'step_thickness': Toggle whether to update the ice thickness (true) or not (false) at each timestep
"""
function TimesteppingParams(;
                        niter0 = 0,
                        dt = 1.0,
                        ntimesteps_velocity_update = 1,
                        end_time = nothing,
                        n_iter_total = nothing, 
                        t0 = nothing,
                        chkpt_freq = Inf,
                        pchkpt_freq = Inf,
                        chkpt_path = "./",
                        step_thickness = true,
                        verbose = false)

    #initialize t0 (really you should read start time from pickup file)
    t0 = niter0 > 0 ? niter0 * dt : 0 
    t0 = map(typeof(dt), t0)

    #check that thickness_timestep_fraction is a positive integer greater than one
    if !(ntimesteps_velocity_update isa Int && ntimesteps_velocity_update >= 1)
        throw(ArgumentError("ntimesteps_velocity_update must be an integer greater than or equal to one"))
    end


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

    return TimesteppingParams(niter0, dt, ntimesteps_velocity_update,end_time, t0, chkpt_freq, pchkpt_freq, 
                            chkpt_path,n_iter_total, n_iter_chkpt, n_iter_pchkpt, step_thickness,verbose)
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