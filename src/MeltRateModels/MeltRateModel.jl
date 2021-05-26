
export AnalyticMeltRate
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


#function get_basal_melt(melt_model::GeometricMeltRateModel, model)
    
 
                            
                            