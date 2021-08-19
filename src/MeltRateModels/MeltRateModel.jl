export AnalyticMeltRate, BinfileMeltRate, add_melt_rate_model!, PlumeEmulator

#add each of the individual melt rate models
include("analytic_melt_rate_model.jl")
include("binfile_melt_rate.jl")
include("plume_emulator.jl")


"""
    function add_melt_rate_model!(model, melt_rate_model)

Interface to endow the model with a melt rate model
"""
function add_melt_rate_model!(model, melt_rate_model)
    ("melt_rate_model" in keys(model.extra_physics)) ? (@info "Model already contained a melt rate model. Overwritten to that just specified...") : model.extra_physics["melt_rate_model"] = melt_rate_model

    #initialize the size of the melt rate grid 
    melt_rate_model = @set melt_rate_model.melt_rate = zeros(size(model.grid.xxh))

    #size(melt_rate_model.melt_rate) == (model.grid.nx, model.grid.ny) ?  model.extra_physics["melt_rate_model"] = melt_rate_model  : @warn "Melt rate grid size does not match model grid size....\n ignoring this melt rate model"
    
    #update the melt rate in both melt model and model
    update_height_above_floatation!(model)
    update_grounded_fraction_on_huv_grids!(model) 
    update_basal_melt!(model)

    return nothing
end
                     

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

