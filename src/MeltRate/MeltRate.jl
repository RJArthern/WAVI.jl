#add each of the individual melt rate models
#include("./analytic_melt_rate_model.jl")
include("./binfile_melt_rate.jl")
include("./plume_emulator.jl")
include("./pico.jl")
             

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


struct MISMIPMeltRateOne{T <: Real} <: AbstractMeltRate 
    α  :: T
    ρi :: T
    ρw :: T
end


"""
    function MISMIPMeltRateOne(; <kwargs>)

Construct a melt rate to specify the melt rate according to the MISMIP_1r experiment

Keyword arguments
=================
- α : Calibrated melt rate prefactor (defaults to unity, value used in MISMIP+)
- ρi : Ice density
- ρw : Water density
"""
MISMIPMeltRateOne(; α = 1.0, ρi = 918.0, ρw = 1028.0) = MISMIPMeltRateOne(α,ρi, ρw)

"""
    update_melt_rate(melt_rate::MISMIPMeltRateOne, fields, grid) 

Update the melt rate when for the MISMIPMeltRateOne type, used in MISMIP Ice 1x experiments
"""
function update_melt_rate!(melt_rate::MISMIPMeltRateOne, fields, grid)
    @unpack basal_melt, h, b  = fields.gh
    draft = -(melt_rate.ρi / melt_rate.ρw) .* h
    cavity_thickness = draft .- b
    cavity_thickness = max.(cavity_thickness, 0)
    m =  melt_rate.α .* 0.2*tanh.(cavity_thickness./75).*max.((-100 .- draft), 0)
    basal_melt[:] .= m[:]
end