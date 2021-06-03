
struct TimesteppingParams{T <: Real, N <: Integer, TO, C, P}
    niter0::N           #initial iteration number
    dt::T                #timestep
    end_time::T          #end time of this simulation
    t0::T                #start time of this simulation 
    chkpt_freq::T        #temporary checkpoint frequency
    pchkpt_freq::T       #permanent checkpoint frequency  
    n_iter_total::TO     #total number of timesteps counting from zero
    n_iter_chkpt::C      #number of iterations per temporary checkpoint
    n_iter_pchkpt::P     #number of iterations per permanent checkpoint
    step_thickness::Bool #toggle whether to step the thickness at each timestep or not
end

function TimesteppingParams(;
                        niter0 = 0,
                        dt = 1.0,
                        end_time = 1.0,
                        t0 = nothing,
                        chkpt_freq = Inf,
                        pchkpt_freq = Inf,
                        step_thickness = true)


    #initialize t0 (really you should read start time from pickup file)
    t0 = niter0 > 0 ? niter0 * dt : 0 
    t0 = map(typeof(dt), t0)

    #compute number of timesteps (total and per checkpoint)
    end_time == Inf ? n_iter_total = Inf : n_iter_total  = round(Int, end_time/dt)
    chkpt_freq == Inf ? n_iter_chkpt = Inf : n_iter_chkpt  = round(Int, chkpt_freq/dt)
    pchkpt_freq == Inf ? n_iter_pchkpt = Inf : n_iter_pchkpt = round(Int, pchkpt_freq/dt)
    
    return TimesteppingParams(niter0, dt, end_time, t0, chkpt_freq, pchkpt_freq,
                                 n_iter_total, n_iter_chkpt, n_iter_pchkpt, step_thickness)
end