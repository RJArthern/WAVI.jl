
#Struct to hold model parameters.
#Format: fieldname::Type = default_value.
#T & N are type parameters, usually real numbers (e.g. Float64) and integers (e.g. Int64) respectively.

#bed_elevation::Array{T,2} = zeros(nx,ny); @assert size(bed_elevation)==(nx,ny)
#starting_thickness::Array{T,2} = zeros(nx,ny); @assert size(starting_thickness)==(nx,ny)

@with_kw struct Params{T <: Real}
    g::T = 9.81
    density_ice::T = 918.0
    density_ocean::T = 1028.0
    gas_const=8.314;
    sec_per_year::T = 3600*24*365.25
    default_thickness::T = 100.
    default_viscosity::T = 1.0e7
    default_temperature::T = 265.700709
    default_damage::T = 0.0
    accumulation_rate::T = 0.0
    basal_melt_rate::T = 0.0
    glen_a_activation_energy::T = 5.8631e+04
    glen_a_ref::T = 4.9e-16 *sec_per_year * 1.0e-9
    glen_temperature_ref::T = 263.15
    glen_n::T = 3.0
    glen_reg_strain_rate::T = 1.0e-5
    weertman_c::T = 1e4
    weertman_m::T = 3.0
    weertman_reg_speed::T = 1.0e-5
    sea_level_wrt_geoid::T = 0.0
    minimum_thickness::T = 50.0
    end
    
    
#structure to hold the solver parameters
@with_kw struct SolverParams{T <: Real, N <: Integer}
    n_iter_viscosity::N = 2;  @assert n_iter_viscosity ==2
    tol_picard::T = 1e-5
    maxiter_picard::N = 30
    tol_coarse::T = 1e-5
    maxiter_coarse::N = 1000
    levels::N = 3
    wavelet_threshold::T = 10.0
    nsmooth::N = 5
    smoother_omega::T = 1.0
    stencil_margin::N = 3
end

struct TimesteppingParams{T <: Real, N <: Integer, TO, C, P}
    n_iter0::N      #initial iteration number
    dt::T           #timestep
    end_time::T     #end time of this simulation
    t0::T           #start time of this simulation 
    chkpt_freq::T   #temporary checkpoint frequency
    pchkpt_freq::T  #permanent checkpoint frequency  
    n_iter_total::TO #total number of timesteps counting from zero
    n_iter_chkpt::C #number of iterations per temporary checkpoint
    n_iter_pchkpt::P#number of iterations per permanent checkpoint
end

function TimesteppingParams(;
                        n_iter0 = 0,
                        dt = 1.0,
                        end_time = 1.0,
                        t0 = nothing,
                        chkpt_freq = Inf,
                        pchkpt_freq = Inf)

    #(n_iter0 > 0) || ArgumentError("n_iter0 must be a positive number")

    #if n_iter0 > 0, check file exists and get start time, else throw error 


    #initialize t0 (really you should read start time from pickup file)
    t0 = n_iter0 > 0 ? n_iter0 * dt : 0 
    t0 = map(typeof(dt), t0)

    #compute number of timesteps (total and per checkpoint)
    end_time == Inf ? n_iter_total = Inf : n_iter_total  = round(Int, end_time/dt)
    chkpt_freq == Inf ? n_iter_chkpt = Inf : n_iter_chkpt  = round(Int, chkpt_freq/dt)
    pchkpt_freq == Inf ? n_iter_pchkpt = Inf : n_iter_pchkpt = round(Int, pchkpt_freq/dt)
    
    return TimesteppingParams(n_iter0, dt, end_time, t0, chkpt_freq, pchkpt_freq, n_iter_total, n_iter_chkpt, n_iter_pchkpt)
end
