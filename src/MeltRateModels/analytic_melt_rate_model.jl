struct AnalyticMeltRate{PC, F, A, M} <: AbstractMeltRateModel{PC,M}
    melt_partial_cell::PC        #specify whether melt applies to partial cells or not
    melt_rate_function::F        #Melt rate function
    function_arguments::A        #Arguments of the function (must be a named tuple, ordered according to function arguments)
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
        melt_rate = zeros(size(melt_rate)) #let melt rate be initialized when we add it to the model
    catch 
        throw(ArgumentError("melt rate function input does not have a method that matches (did you remember to vectorize your melt_rate_function?)"))
    end
    
    return AnalyticMeltRate(melt_partial_cell, melt_rate_function, function_arguments, melt_rate)
end

"""
    function update_melt_rate_model!(melt_model::AnalyticMeltRate, model)

Update the melt rate for an analytic-type melt rate model
"""
function update_melt_rate_model!(melt_model::AnalyticMeltRate, model)
    @unpack melt_rate, function_arguments, melt_rate_function = melt_model
        if melt_model.melt_partial_cell
        melt_rate .= melt_rate_function(function_arguments...).* (1 .- model.fields.gh.grounded_fraction)
    else
        grounded_matrix =  model.fields.gh.grounded_fraction .- mod.(model.fields.gh.grounded_fraction,1) #send any partially floating cells to zero grounded fraction
        melt_rate .= melt_rate_function(function_arguments...).* (1 .- grounded_matrix)
    end
    return nothing
end
    
