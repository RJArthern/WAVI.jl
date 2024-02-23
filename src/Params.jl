struct Params{T <: Real, A, W, G}
                      dt :: T
                       g :: T
             density_ice :: T
           density_ocean :: T
               gas_const :: T 
           sec_per_year  :: T
       default_thickness :: T 
       default_viscosity :: T 
     default_temperature :: T
          default_damage :: T 
       accumulation_rate :: A
glen_a_activation_energy :: T 
              glen_a_ref :: G
    glen_temperature_ref :: T 
                  glen_n :: T 
    glen_reg_strain_rate :: T 
              weertman_c :: W
              weertman_m :: T 
      weertman_reg_speed :: T 
     sea_level_wrt_geoid :: T
       minimum_thickness :: T 
                step_haf :: Bool   
           evolveShelves :: Bool
                smallHAF :: T
end


"""
Params(; <kwargs>)

Construct a WAVI.jl parameters object for holding physical parameters.

Keyword arguments
=================
- `dt`: model timestep (NB: simulation timestep set in timestepping params, this value is updated when model embedded to the value specified in timestepping_params when passed to simulation)
- `g`: gravitational acceleration (m^2 / s)
- `density_ice`: ice density (kg / m^3)
- `density_ocean`: ocean water density (kg / m^3)
- `gas_const`: gas constant in glen b calculation
- `sec_per_year`: seconds per year (s)
- `default_thickness`: thickness value reverted to if no initial thickness passed (m)
- `default_viscosity`: viscosity value reverted to if no initial thickness passed (Pa s)
- `default_temperature`: temperature value reverted to if no initial thickness passed (K)
- `default_damage`: damage value reverted to if no initial thickness passed (dimensionless)
- `accumulation_rate`: uniform accumulation_rate (m/yr)
- `glen_a_activation_energy`: activation energy in glen b calculation
- `glen_a_ref`: array of glen a reference values used in glen b calculation
- `glen_temperature_ref`: reference temperature using in glen b calculation
- `glen_n`: exponent in glen b calculation
- `glen_reg_strain_rate`: strain rate regularization value
- `weertman_c`: basal sliding field of coefficients
- `weertman_m`: sliding law exponent
- `weertman_reg_speed`: regularization speed, used to prevent bed speed going to zero
- `sea_level_wrt_geoid`: reference sea level
- `minimum_thickness`: minimum ice thickness on model domain
- `evolveShelves`: flag for turning on and off the evolution of the shelves in the forward run_simulation
- `smallHAF`: small value of HAF used within update_thickness when not evolving shelves
"""
function Params(; g = 9.81, 
                  density_ice = 918.0,
                  density_ocean = 1028.0,
                  gas_const = 8.314, 
                  sec_per_year =3600*24*365.25,
                  default_thickness= 100.,
                  default_viscosity= 1.0e7,
                  default_temperature=263.15,
                  default_damage= 0.0,
                  accumulation_rate= 0.0,
                  glen_a_activation_energy = 5.8631e+04,
                  glen_a_ref= 4.9e-16 *sec_per_year * 1.0e-9,
                  glen_temperature_ref= 263.15,
                  glen_n = 3.0,
                  glen_reg_strain_rate = 1.0e-5,
                  weertman_c = 1.0e4,
                  weertman_m  = 3.0,
                  weertman_reg_speed = 1.0e-5,
                  sea_level_wrt_geoid  = 0.0,
                  minimum_thickness = 50.0,
                  step_haf = true,
                  evolveShelves = true,
                  smallHAF = 1.0)
                      
  #defualt the timestep to 1.0 (will be updated when the model is embedded in a simulation)
  dt = 1.0

  return Params(
                  dt, 
                  g, 
                  density_ice, 
                  density_ocean, 
                  gas_const,
                  sec_per_year, 
                  default_thickness, 
                  default_viscosity,
                  default_temperature,
                  default_damage,
                  accumulation_rate,
                  glen_a_activation_energy,
                  glen_a_ref,
                  glen_temperature_ref,
                  glen_n,
                  glen_reg_strain_rate,
                  weertman_c,
                  weertman_m,
                  weertman_reg_speed,
                  sea_level_wrt_geoid,
                  minimum_thickness,
                  step_haf,
                  evolveShelves,
                  smallHAF
                  )
end
