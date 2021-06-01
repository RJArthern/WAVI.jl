#structure to hold the solver parameters
@with_kw struct Params{T <: Real}
dt::T = 1.0
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