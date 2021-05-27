
export AnalyticMeltRate, add_melt_rate_model!
#include("./")

struct AnalyticMeltRate{PC, F, A, M} <: AbstractMeltRateModel{PC,M}
    melt_partial_cell::PC        #specify whether melt applies to partial cells or not
    melt_rate_function::F        #Melt rate function
    function_arguments::A        #Arguments of the function (named tuple)
    melt_rate::M                 #Current melt rate
end

function AnalyticMeltRate(;  
                            melt_partial_cell = false,
                            melt_rate_function = nothing,
                            function_arguments = nothing)

    #check that a melt rate function has been passed
    ~(melt_rate_function === nothing) || throw(ArgumentError("You must pass an array to this function..."))
    
    #check that function works, and, if so, assign to melt rate
    melt_rate = nothing
    try melt_rate_function(function_arguments...) 
        melt_rate = melt_rate_function(function_arguments...) 
    catch 
        throw(ArgumentError("melt rate function input does not have a method that matches"))
    end
    
    return AnalyticMeltRate(melt_partial_cell, melt_rate_function, function_arguments, melt_rate)
end


function update_melt_rate_model!(melt_model::AnalyticMeltRate, model)
    @unpack melt_rate, function_arguments, melt_rate_function = melt_model
    melt_rate .= melt_rate_function(function_arguments...) #this is not actually necessary because the melt rate function updates automatically when arguments are updated.. 
    return nothing
end
    


"""
    function add_melt_rate_model!(model, melt_rate_model)

Interface to endow the model with a melt rate model
"""
function add_melt_rate_model!(model, melt_rate_model)
    ~("melt_rate_model" in keys(model.extra_physics)) || @info "Model already contained a melt rate model. Overwritten to that just specified..."
    size(melt_rate_model.melt_rate) == (model.grid.nx, model.grid.ny) ?  model.extra_physics["melt_rate_model"] = melt_rate_model  : @warn "Melt rate grid size does not match model grid size....\n ignoring this melt rate model"
   
    return nothing
end
                            
                            