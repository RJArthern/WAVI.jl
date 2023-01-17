struct MeltRateExponentVariationBasins{T <: Real} <: AbstractMeltRate
    γT1 :: T                       #(calibrated) heat exchange velocity in basin 1
    γT2 :: T                       #(calibrated) heat exchange velocity in basin 2
    basinID_1 :: T                 # basin ID value for basin 1
    basinID_2 :: T                 # basin ID value for basin 2
    λ1 :: T                       #liquidus slope 
    λ2 :: T                       #liquidus intercept 
    λ3 :: T                       #liquidus pressure coefficient
    ρi :: T                       #ice density 
    ρw :: T                       #water density
    L :: T                        #latent heat of fusion for ice 
    c :: T                        #specific heat capacity of ocean water 
    Ta :: Function                #Function specifying ambient temperature  
    Sa :: Function                #Function specifying ambient salinity
    flocal :: Bool                #Flag for local or non-local temperature dependence
    melt_partial_cell :: Bool       #Flag for melting applied to partial cells or not
    melt_exp :: T                  # exponent in the melt rate function
end

"""
MeltRateExponentVariationBasins(; ,kwargs)

Construct a MeltRateExponentVariationBasins object to prescribe the melt rate in WAVI when the domain is split into two basins

Keyword arguments 
=================
- γT1: calibrated heat exchange velocity in basin 1 (units m/s).
- γT2: calibrated heat exchange velocity in basin 2 (units m/s).
- basinID_1: ID number of basin 1 
- basinID_2: ID number of basin 2 
- λ1: liquidus slope (units ∘C)
- λ2: liquidus intercept (units ∘C)
- λ3: liquidus pressure coefficient (units K/m)
- ρi: ice density (units kg/m^3)
- ρw: water density (units kg/m^3)
- L: latent heat of fusion for ice (units: J/kg)
- c: specific heat capacity of water (units J/kg/K)
- Ta: ambient temperature function, defaults to the ISOMIP warm0 scenario. Ta must be a function of a single variable (depth) [time dependence not yet implemented in WAVI.jl]
- Sa: ambient salnity function, defaults to the ISOMIP warm0 scenario.  Sa must be a function of a single variable (depth) [time dependence not yet implemented in WAVI.jl]
- flocal: flag for local or non-local temperature dependence. In this version of the function, only local dependence is allowed  (eqn 4 in Favier 2019 GMD) 
- melt_partial_cell: flag for specify whether to apply melt on partially floating cells (true) or not (false)
- melt_exp: exponent in the melt rate function
"""

function MeltRateExponentVariationBasins(;
                            γT1 = 1.e-3,
                            γT2 = 1.e-3,
                            basinID_1=1.0,
                            basinID_2=2.0,
                            λ1 = -5.73e-2,
                            λ2 = 8.32e-4,
                            λ3 = 7.61e-4,
                            L = 3.35e5,
                            c = 3.974e3,
                            ρi = 918.8,
                            ρw = 1028.0,
                            Ta = isomip_warm0_temp,
                            Sa = isomip_warm0_salinity,
                            flocal = true,
                            melt_partial_cell = false,
                            melt_exp = 2.0)

    #check that Ta and Sa accept a single argument
    try Ta(0.0)
    catch
        ArgumentError("Ambient temperature function Ta must accept a single argument")
    end
    try Sa(0.0)
    catch
        ArgumentError("Ambient salinity function Sa must accept a single argument")
    end

    #this function is only value for flocal = true:
    @assert flocal "This function is only valid for a local temperature dependece (flocal=true)"

    return MeltRateExponentVariationBasins(γT1, γT2, basinID_1, basinID_2, λ1, λ2, λ3, ρi, ρw, L, c, Ta, Sa, flocal, melt_partial_cell, melt_exp)
end

"""
    update_melt_rate!(quad_melt_rate::MeltRateExponentVariationBasins, fields, grid)

Wrapper script to update the melt rate for a MeltRateExponentVariationBasins.
"""
function update_melt_rate!(quad_melt_rate::MeltRateExponentVariationBasins, fields, grid, clock)
    @unpack basal_melt, h, b, grounded_fraction = fields.gh #get the ice thickness and grounded fraction
 
    #compute the ice draft
    zb = fields.gh.b .* (grounded_fraction .== 1) + - quad_melt_rate.ρi / quad_melt_rate.ρw .* fields.gh.h .* (grounded_fraction .< 1)

    basinIDs=grid.basin_ID   
    set_melt_rate_exponent_variation_basins!(basal_melt,
                            quad_melt_rate,
                            zb, 
                            grounded_fraction,
                            basinIDs)
    return nothing
end


function set_melt_rate_exponent_variation_basins!(basal_melt,
                                qmr,
                                zb, 
                                grounded_fraction,
                                basinIDs)

    #compute freezing temperature and Ta - Tf everywhere
    Sa_shelf = qmr.Sa.(zb)
    Ta_shelf = qmr.Ta.(zb)
#    Tstar = qmr.λ1 .* Sa_shelf .+ qmr.λ2 .+ qmr.λ3 .* zb .- Ta_shelf
    Tstar =  Ta_shelf .- (qmr.λ1 .* Sa_shelf .+ qmr.λ2 .+ qmr.λ3 .* zb)
    Tstar_shelf_mean = sum(Tstar[grounded_fraction .== 0])/length(Tstar[grounded_fraction .== 0])

    #basal_melt = zeros(size(grounded_fraction));

    #set melt rate
    if (qmr.melt_partial_cell) && (qmr.flocal) #partial cell melting and local 
        basal_melt[basinIDs .==qmr.basinID_1].=  qmr.γT1 .* (qmr.ρw * qmr.c / qmr.ρi /qmr.L)^(qmr.melt_exp) .* Tstar[basinIDs .==qmr.basinID_1].^(qmr.melt_exp) 
        basal_melt[basinIDs .==qmr.basinID_2].=  qmr.γT2 .* (qmr.ρw * qmr.c / qmr.ρi /qmr.L)^(qmr.melt_exp) .* Tstar[basinIDs .==qmr.basinID_2].^(qmr.melt_exp) 
        basal_melt[:] .=  basal_melt[:] .* (1 .- grounded_fraction[:])

    elseif ~(qmr.melt_partial_cell) && (qmr.flocal) #no partial cell melting and local
        basal_melt[basinIDs .==qmr.basinID_1] .=  qmr.γT1 .* (qmr.ρw * qmr.c / qmr.ρi /qmr.L)^(qmr.melt_exp).* Tstar[basinIDs .==qmr.basinID_1].^(qmr.melt_exp)
        basal_melt[basinIDs .==qmr.basinID_2] .=  qmr.γT2 .* (qmr.ρw * qmr.c / qmr.ρi /qmr.L)^(qmr.melt_exp).* Tstar[basinIDs .==qmr.basinID_2].^(qmr.melt_exp)
        basal_melt[.~(grounded_fraction .== 0)] .= 0 
    end
    basal_melt[:] .= basal_melt[:].* 365.25*24*60*60
    return nothing
end
