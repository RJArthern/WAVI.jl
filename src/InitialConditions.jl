"""
InitialConditions(; 
                    initial_thickness = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_viscosity = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_temperature = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_damage = fill!(Array{Float64}(undef,1,1),NaN)),
                    initial_u = fill!(Array{Float64}(undef,1,1),NaN)),
                    initial_v = fill!(Array{Float64}(undef,1,1),NaN))

Construct a WAVI.jl initial conditions object. 
Unpassed arguments default to 1x1 nan matrix; unspecified initial conditions are overwritten by default values specified in Params structure when model is assembled.

Keyword arguments
=================

- 'initial_thickness': (nx x ny) matrix defining ice thickness at t = 0
- 'initial_viscosity': (nx x ny x nz) matrix defining viscosity on sigma levels at t = 0
- 'initial_temperature': (nx x ny) matrix defining temperature on sigma levels at t = 0
- 'initial_damage': (nx x ny) matrix defining depth averaged damage at t = 0
- 'initial_u': (nxu x nyu) matrix defining ice velocity in x direction at t = 0
- 'initial_v': (nxv x nyv) matrix defining ice velocity in y direction at t = 0
"""
@with_kw struct InitialConditions{T <: Real}
    initial_thickness::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_viscosity::Array{T,3} = fill!(Array{Float64}(undef,1,1,1),NaN)
    initial_temperature::Array{T,3} = fill!(Array{Float64}(undef,1,1,1),NaN)
    initial_damage::Array{T,3} = fill!(Array{Float64}(undef,1,1,1),NaN)
    initial_u::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_v::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
end
