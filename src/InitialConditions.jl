#structure to store initial conditions
@with_kw struct InitialConditions{T <: Real}
    initial_thickness::Array{T,2} = zeros(10,10)
    initial_viscosity::Array{T,2} = zeros(10,10) #placeholder array
    initial_temperature::Array{T,2} = zeros(10,10)
    initial_damage::Array{T,2} = zeros(10,10)
end