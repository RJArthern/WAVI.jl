#structure to store initial conditions. We default values are a 1x1 nan matrix. Unspecified initial conditions are overwritten in Model constructor to default values passed via Params structure
@with_kw struct InitialConditions{T <: Real}
    initial_thickness::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_viscosity::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_temperature::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_damage::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
end