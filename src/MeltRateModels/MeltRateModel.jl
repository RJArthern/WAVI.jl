export AnalyticMeltRate, BinfileMeltRate, add_melt_rate_model!

#add each of the individual melt rate models
include("analytic_melt_rate_model.jl")
include("binfile_melt_rate.jl")

"""
    function add_melt_rate_model!(model, melt_rate_model)

Interface to endow the model with a melt rate model
"""
function add_melt_rate_model!(model, melt_rate_model)
    ~("melt_rate_model" in keys(model.extra_physics)) || @info "Model already contained a melt rate model. Overwritten to that just specified..."
    size(melt_rate_model.melt_rate) == (model.grid.nx, model.grid.ny) ?  model.extra_physics["melt_rate_model"] = melt_rate_model  : @warn "Melt rate grid size does not match model grid size....\n ignoring this melt rate model"
    
    #update the melt rate in both melt model and model
    update_grounded_fraction_on_huv_grids!(model) 
    update_basal_melt!(model)

    return nothing
end
                            
                            