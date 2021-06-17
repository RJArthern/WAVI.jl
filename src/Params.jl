"""
Params(; 
        dt = 1.0,
        g  = 9.81,
        density_ice = 918.0,
        density_ocean = 1028.0,
        gas_const = 8.314, 
        sec_per_year = 3.15576e7,
        default_thickness = 100.0,
        default_viscosity = 1.0e7,
        default_temperature = 265.700709,
        default_damage = 0.0,
        accumulation_rate = 0.0,
        glen_a_activation_energy = 5.8631e+04,
        glen_a_ref = 4.9e-16 *sec_per_year * 1.0e-9,
        glen_temperature_ref = 263.15,
        glen_n = 3.0,
        glen_reg_strain_rate = 1.0e-5,
        weertman_c = 1.0e4
        weertman_m = 3.0,
        weertman_reg_speed = 1.0e-5, 
        sea_level_wrt_geoid = 0.0,
        minimum_thickness = 50.0)

Construct a WAVI.jl parameters object for holding physical parameters.

Keyword arguments
=================
- 'dt': model timestep (NB: simulation timestep set in timestepping params, this value is updated when model embedded to the value specified in timestepping_params when passed to simulation)
- 'g': gravitational acceleration (m^2 / s)
- 'density_ice': ice density (kg / m^3)
- 'density_ocean' ocean water density (kg / m^3)
- 'gas_const': gas constant in glen b calculation
- 'sec_per_year': seconds per year (s)
- 'default_thickness': thickness value reverted to if no initial thickness passed (m)
- 'default_viscosity': viscosity value reverted to if no initial thickness passed (Pa s)
- 'default_temperature': temperature value reverted to if no initial thickness passed (K)
- 'default_damage': damage value reverted to if no initial thickness passed (dimensionless)
- 'accumulation_rate': uniform accumulation_rate (m/yr)
- 'glen_a_activation_energy': activation energy in glen b calculation
- 'glen_a_ref': glen a reference value used in glen b calculation
- 'glen_temperature_ref': reference temperature using in glen b calculation
- 'glen_n': exponent in glen b calculation
- 'glen_reg_strain_rate': strain rate regularization value
- 'weertman_c': uniform basal sliding coefficient
- 'weertman_m': sliding law exponent
- 'weertman_reg_speed': regularization speed, used to prevent bed speed going to zero
- 'sea_level_wrt_geoid': reference sea level
- 'minimum_thickness': minimum ice thickness on model domain
"""
@with_kw struct Params{T <: Real}
                      dt :: T = 1.0
                       g :: T = 9.81
             density_ice :: T = 918.0
           density_ocean :: T = 1028.0
               gas_const :: T = 8.314;
           sec_per_year  :: T = 3600*24*365.25
       default_thickness :: T = 100.
       default_viscosity :: T = 1.0e7
     default_temperature :: T = 265.700709
          default_damage :: T = 0.0
       accumulation_rate :: T = 0.0
glen_a_activation_energy :: T = 5.8631e+04
              glen_a_ref :: T = 4.9e-16 *sec_per_year * 1.0e-9
    glen_temperature_ref :: T = 263.15
                  glen_n :: T = 3.0
    glen_reg_strain_rate :: T = 1.0e-5
              weertman_c :: T = 1e4
              weertman_m :: T = 3.0
      weertman_reg_speed :: T = 1.0e-5
     sea_level_wrt_geoid :: T = 0.0
       minimum_thickness :: T = 50.0
end