#add each of the individual melt rate models
#include("./analytic_melt_rate_model.jl")
include("./binfile_melt_rate.jl")
include("./plume_emulator.jl")
include("./pico.jl")
include("./quadratic_melt_rate.jl")
include("./quadratic_melt_rate_exponent_variation.jl")
include("./mismip_melt_rate.jl")
include("./uniform_melt_float_only.jl")
include("./uniform_melt_float_only_basin_specific.jl")

             

##### default temperature and salinity profiles #####
"""
    two_layer_ambient()

Returns a function that, in turn, returns a value v_low below a depth d_low, and a value v_hi above a depth d_hi > d_low, with a linear interpolation between these value.
"""
function two_layer_function(z;v_low,v_hi,d_low,d_hi)
    @assert d_low < d_hi
    if  z< d_low
        return v_low
    elseif z > d_hi
        return v_hi
    else #depth between d_low and d_hi
        return v_low + (v_hi - v_low)/(d_hi - d_low) * (z - d_low)
    end
end


isomip_warm0_salinity(z) = two_layer_function(z, v_low = 34.6, v_hi = 34.0, d_low = -700, d_hi = -300)
isomip_warm0_temp(z) = two_layer_function(z, v_low =1.2, v_hi = -1.0, d_low = -700, d_hi = -300)


struct UniformMeltRate{T <: Real} <: AbstractMeltRate 
    m :: T #uniform melt rate applied everywhere
end

UniformMeltRate(; m = 0.0) = UniformMeltRate(m) 

"""
    update_melt_rate(melt_rate::UniformMeltRate, fields, grid) 

Update the melt rate when for the UniformMeltRate type

"""
function update_melt_rate!(melt_rate::UniformMeltRate, fields, grid) 
    @unpack basal_melt = fields.gh
    basal_melt .= melt_rate.m
    return nothing
end

