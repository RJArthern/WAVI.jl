module Parameters

import WAVI.Grids: reconstruct_on_grid, reconstruct_on_subdomain
export reconstruct_on_grid, reconstruct_on_subdomain
export Params, SolverParams, TimesteppingParams

using Parameters

using WAVI.Grids
using WAVI.Time

struct Params{T, A <: Union{T,Array{T,2}}, G <: Union{T,Array{T,2}}, H <: Union{T,Array{T,2}}, E <: Union{T,Array{T,2}}}
                      dt :: T
                       g :: T
             density_ice :: T
      density_freshwater :: T
           density_ocean :: T
               gas_const :: T 
           sec_per_year  :: T
       default_thickness :: T 
       default_viscosity :: T 
     default_temperature :: T
          default_damage :: T 
  default_strain_history :: T
    ice_tensile_strength :: T
ice_compressive_strength :: T
 critical_elastic_energy :: T
      phase_field_length :: T
     energy_release_rate :: T
degradation_regularisation :: T
       accumulation_rate :: A
glen_a_activation_energy :: T 
              glen_a_ref :: G
    glen_temperature_ref :: T 
                  glen_n :: T 
    glen_reg_strain_rate :: T 
          elastic_lambda :: T 
              elastic_mu :: T 
     sea_level_wrt_geoid :: T
       minimum_thickness :: T 
           evolveShelves :: Bool
                smallHAF :: T
   basal_water_thickness :: T
   hydraulic_potential_b :: H
      effective_pressure :: E
 default_temperature_ave :: T
      default_preBfactor :: T
end


"""
Params(; <kwargs>)

Construct a WAVI.jl parameters object for holding physical parameters.

Keyword arguments
=================
- `dt`: model timestep (NB: simulation timestep set in timestepping params, this value is updated when model embedded to the value specified in timestepping_params when passed to simulation)
- `g`: gravitational acceleration (m / s^2)
- `density_ice`: ice density (kg / m^3)
- `density_freshwater`: freshwater density (kg / m^3)
- `density_ocean`: ocean water density (kg / m^3)
- `gas_const`: gas constant in glen b calculation
- `sec_per_year`: seconds per year (s)
- `default_thickness`: thickness value reverted to if no initial thickness passed (m)
- `default_viscosity`: viscosity value reverted to if no initial thickness passed (Pa s)
- `default_temperature`: temperature value reverted to if no initial thickness passed (K)
- `default_damage`: damage value reverted to if no initial damage passed (dimensionless)
- `default_strain_history`: maximum previous strain energy reverted to if no initial value passed (Pa)
- `ice_tensile_strength`: uniaxial tensile strength of ice (Pa)
- `ice_compressive_strength`: uniaxial conpressive strength of ice (Pa)
- `critical_elastic_energy`: critical elastic strain energy for fracture nucleation (Pa)
- `phase_field_length`: phase field length for fracture (m)
- `energy_release_rate`: energy release rate for fracture (N/m)
- `degradation_regularisation`: smallest degradation allowed
- `accumulation_rate`: uniform accumulation_rate (m/yr)
- `glen_a_activation_energy`: activation energy in glen b calculation
- `glen_a_ref`: array of glen a reference values used in glen b calculation
- `glen_temperature_ref`: reference temperature using in glen b calculation
- `glen_n`: exponent in glen b calculation
- `glen_reg_strain_rate`: strain rate regularization value
- `elastic_lambda`: First elastic Lame parameter (Pa)
- `elastic_mu`: Second elastic Lame parameter (Pa) 
- `sea_level_wrt_geoid`: reference sea level
- `minimum_thickness`: minimum ice thickness on model domain
- `evolveShelves`: flag for turning on and off the evolution of the shelves in the forward run_simulation
- `smallHAF`: small value of HAF used within update_thickness when not evolving shelves
- `basal_water_thickness` : basal water thickness (m)
- `hydraulic_potential_b` : hydraulic_potential at the bed (Pa)
- `effective_pressure` : effective pressure (Pa)
- `default_temperature_ave` : depth-averaged temperature (K)
- `default_preBfactor` : default viscosity enhancement factor (non-dimensional).
"""
function Params(; g = 9.81, 
                  density_ice = 918.0,
                  density_freshwater = 1000.0,
                  density_ocean = 1028.0,
                  gas_const = 8.314, 
                  sec_per_year =3600*24*365.25,
                  default_thickness= 100.,
                  default_viscosity= 1.0e7,
                  default_temperature=263.15,
                  default_damage= 0.0,
                  default_strain_history = 0.0,
                  ice_tensile_strength = 262.9e3,
                  ice_compressive_strength =  496.7e3,
                  critical_elastic_energy = 1.0,
                  phase_field_length = 1.0,
                  energy_release_rate = 1.0,
                  degradation_regularisation = 1e-2,
                  accumulation_rate= 0.0,
                  glen_a_activation_energy = 5.8631e+04,
                  glen_a_ref= 4.9e-16 *sec_per_year * 1.0e-9,
                  glen_temperature_ref= 263.15,
                  glen_n = 3.0,
                  glen_reg_strain_rate = 1.0e-5,
                  elastic_lambda = 6.54e9,
                  elastic_mu = 3.52e9,
                  sea_level_wrt_geoid  = 0.0,
                  minimum_thickness = 50.0,
                  evolveShelves = true,
                  smallHAF = 1.0,
                  basal_water_thickness = 0.0,
                  hydraulic_potential_b = 0.0,
                  effective_pressure = 1.0e6,
                  default_temperature_ave = 253.15,
                  default_preBfactor = 1.0)
                      
  #default the timestep to 1.0 (will be updated when the model is embedded in a simulation)
  dt = 1.0

  return Params(
                  dt, 
                  g, 
                  density_ice,
                  density_freshwater, 
                  density_ocean, 
                  gas_const,
                  sec_per_year, 
                  default_thickness, 
                  default_viscosity,
                  default_temperature,
                  default_damage,
                  default_strain_history,
                  ice_tensile_strength,
                  ice_compressive_strength,
                  critical_elastic_energy,
                  phase_field_length,
                  energy_release_rate,
                  degradation_regularisation,
                  accumulation_rate,
                  glen_a_activation_energy,
                  glen_a_ref,
                  glen_temperature_ref,
                  glen_n,
                  glen_reg_strain_rate,
                  elastic_lambda,
                  elastic_mu,
                  sea_level_wrt_geoid,
                  minimum_thickness,
                  evolveShelves,
                  smallHAF,
                  basal_water_thickness,
                  hydraulic_potential_b,
                  effective_pressure,
                  default_temperature_ave,
                  default_preBfactor
                  )
end

#Outer constructor that expands scalars into arrays of approriate size.
function reconstruct_on_grid(params::Params, grid::Grid) 
     return     Params(
                  params.dt, 
                  params.g, 
                  params.density_ice,
                  params.density_freshwater, 
                  params.density_ocean, 
                  params.gas_const,
                  params.sec_per_year, 
                  params.default_thickness, 
                  params.default_viscosity,
                  params.default_temperature,
                  params.default_damage,
                  params.default_strain_history,
                  params.ice_tensile_strength,
                  params.ice_compressive_strength,
                  params.critical_elastic_energy,
                  params.phase_field_length,
                  params.energy_release_rate,
                  params.degradation_regularisation,
                  field_on_grid(params.accumulation_rate, grid),
                  params.glen_a_activation_energy,
                  field_on_grid(params.glen_a_ref, grid),
                  params.glen_temperature_ref,
                  params.glen_n,
                  params.glen_reg_strain_rate,
                  params.elastic_lambda,
                  params.elastic_mu,
                  params.sea_level_wrt_geoid,
                  params.minimum_thickness,
                  params.evolveShelves,
                  params.smallHAF,
                  params.basal_water_thickness,
                  field_on_grid(params.hydraulic_potential_b, grid),
                  field_on_grid(params.effective_pressure, grid),
                  params.default_temperature_ave,
                  params.default_preBfactor
                  )
end

#Outer constructor that selects parameters for a specified subdomain from parameters defined on a particular grid
function reconstruct_on_subdomain(params::Params, grid::Grid, subdomain::NTuple{4,<: Integer})
    return Params(  params.dt, 
                    params.g, 
                    params.density_ice,
                    params.density_freshwater, 
                    params.density_ocean, 
                    params.gas_const,
                    params.sec_per_year, 
                    params.default_thickness, 
                    params.default_viscosity,
                    params.default_temperature,
                    params.default_damage,
                    params.default_strain_history,
                    params.ice_tensile_strength,
                    params.ice_compressive_strength,
                    params.critical_elastic_energy,
                    params.phase_field_length,
                    params.energy_release_rate,
                    params.degradation_regularisation,
                    spatial_on_subdomain(params.accumulation_rate, grid, subdomain),
                    params.glen_a_activation_energy,
                    spatial_on_subdomain(params.glen_a_ref, grid, subdomain),
                    params.glen_temperature_ref,
                    params.glen_n,
                    params.glen_reg_strain_rate,
                    params.elastic_lambda,
                    params.elastic_mu,
                    params.sea_level_wrt_geoid,
                    params.minimum_thickness,
                    params.evolveShelves,
                    params.smallHAF,
                    params.basal_water_thickness,
                    spatial_on_subdomain(params.hydraulic_potential_b, grid, subdomain),
                    spatial_on_subdomain(params.effective_pressure, grid, subdomain),
                    params.default_temperature_ave,
                    params.default_preBfactor
                    )
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
    super_implicitness::T = 1.0
end

struct TimesteppingParams{T <: Real, N <: Integer, TO, C, P, Y <: Real}
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
ntimesteps_climate_forcing_update :: N   #number of timesteps per update of the climate forcing
                    ref_time :: Y        #set the reference time of the simulation (e.g. 2015)
end


"""
TimesteppingParams(;
                    niter0 = 0,
                    dt = 1.0,
                    ntimesteps_velocity_update = 1,
                    end_time = 1.0,
                    n_iter_total = nothing, 
                    chkpt_freq = Inf,
                    pchkpt_freq = Inf,
                    chkpt_path = './',
                    step_thickness = true,
                    verbose = false, 
                    ref_time = 0)

Construct a WAVI.jl TimesteppingParams object.
TimesteppingParams stores information relating to timestepping.

Keyword arguments
=================
- 'niter0': Iteration number of the first timestep. niter0 = 0 corresponds to a new simulation, while niter0 > 0 (positive integer) corresponds to a pickup.
- 'dt': Model timestep
- 'ntimesteps_velocity_update': number of substeps at which to update the velocity (i.e. the velocity is updated every dt*ntimesteps_velocity_update)
- 'end_time': Simulation termination time
- 'n_iter_total': Total number of timesteps counting from zero
- 'chkpt_freq': Frequency of outputting temporary checkpoints
- 'pchkpt_freq': Frequecy with which permanent checkpoints are pass
- 'chkpt_path' : Path to location checkpoint output
- 'step_thickness': Toggle whether to update the ice thickness (true) or not (false) at each timestep
- 'verbose': Toggle whether to output when the timestepping have been performed (true) or not (false)
- 'ntimesteps_climate_forcing_update': set the number of timesteps per update of the climate forcing
- 'ref_time': set the reference (clock) time
"""
function TimesteppingParams(;
                        niter0 = 0,
                        dt = 1.0,
                        ntimesteps_velocity_update = 1,
                        end_time = nothing,
                        n_iter_total = nothing, 
                        chkpt_freq = Inf,
                        pchkpt_freq = Inf,
                        chkpt_path = "./",
                        step_thickness = true,
                        verbose = false, 
                        ntimesteps_climate_forcing_update = typemax(Int), 
                        ref_time = 0)

    #initialize t0 (really you should read start time from pickup file)
    t0 = niter0 > 0 ? niter0 * dt : 0 
    t0 = map(typeof(dt), t0)

    #check that thickness_timestep_fraction is a positive integer greater than one
    if !(ntimesteps_velocity_update isa Int && ntimesteps_velocity_update >= 1)
        throw(ArgumentError("ntimesteps_velocity_update must be an integer greater than or equal to one"))
    end


    #check compatibility of n_iter_total and end_time, and compute them 
    end_time, n_iter_total = compute_iterations_and_end_time(end_time, n_iter_total, dt)

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
                            chkpt_path,n_iter_total, n_iter_chkpt, n_iter_pchkpt, step_thickness,verbose, ntimesteps_climate_forcing_update, ref_time)
end

end
