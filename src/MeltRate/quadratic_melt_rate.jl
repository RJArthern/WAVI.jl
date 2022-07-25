struct QuadraticMeltRate{T <: Real} <: AbstractMeltRate
    γT :: T                       #(calibrated) heat exchange velocity
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
end

"""
    QuadraticMeltRate(; ,kwargs)

Construct a QuadraticMeltRate object to prescribe the melt rate in WAVI

Keyword arguments 
=================
- γT: calibrated heat exchange velocity (units m/s).
- λ1: liquidus slope (units ∘C)
- λ2: liquidus intercept (units ∘C)
- λ3: liquidus pressure coefficient (units K/m)
- ρi: ice density (units kg/m^3)
- ρw: water density (units kg/m^3)
- L: latent heat of fusion for ice (units: J/kg)
- c: specific heat capacity of water (units J/kg/K)
- Ta: ambient temperature function, defaults to the ISOMIP warm0 scenario. Ta must be a function of a single variable (depth) [time dependence not yet implemented in WAVI.jl]
- Sa: ambient salnity function, defaults to the ISOMIP warm0 scenario.  Sa must be a function of a single variable (depth) [time dependence not yet implemented in WAVI.jl]
- flocal: flag for local or non-local temperature dependence. Set to true for local dependence (eqn 4 in Favier 2019 GMD) or false for non-local dependence (eqn 5)
- melt_partial_cell: flag for specify whether to apply melt on partially floating cells (true) or not (false)
"""

function QuadraticMeltRate(;
                            γT = 1.e-3,
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
                            melt_partial_cell = false)

    #check that Ta and Sa accept a single argument
    try Ta(0.0)
    catch
        ArgumentError("Ambient temperature function Ta must accept a single argument")
    end
    try Sa(0.0)
    catch
        ArgumentError("Ambient salinity function Sa must accept a single argument")
    end

    return QuadraticMeltRate(γT, λ1, λ2, λ3, ρi, ρw, L, c, Ta, Sa, flocal, melt_partial_cell)
end

"""
    update_melt_rate!(quad_melt_rate::QuadraticMeltRate, fields, grid)

Wrapper script to update the melt rate for a QuadraticMeltRate.
"""
function update_melt_rate!(quad_melt_rate::QuadraticMeltRate, fields, grid)
    @unpack basal_melt, h, b, grounded_fraction = fields.gh #get the ice thickness and grounded fraction
 
    #compute the ice draft
    zb = fields.gh.b .* (grounded_fraction .== 1) + - quad_melt_rate.ρi / quad_melt_rate.ρw .* fields.gh.h .* (grounded_fraction .< 1)

    set_quadratic_melt_rate!(basal_melt,
                            quad_melt_rate,
                            zb, 
                            grounded_fraction)
    return nothing
end


function set_quadratic_melt_rate!(basal_melt,
                                qmr,
                                zb, 
                                grounded_fraction)

    #compute freezing temperature and Ta - Tf everywhere
    Sa_shelf = qmr.Sa.(zb)
    Ta_shelf = qmr.Ta.(zb)
    Tstar = qmr.λ1 .* Sa_shelf .+ qmr.λ2 .+ qmr.λ3 .* zb .- Ta_shelf
    Tstar_shelf_mean = sum(Tstar[grounded_fraction .== 0])/length(Tstar[grounded_fraction .== 0])

    #set melt rate
    if (qmr.melt_partial_cell) && (qmr.flocal) #partial cell melting and local 
        basal_melt[:] .=  qmr.γT .* (qmr.ρw * qmr.c / qmr.ρi /qmr.L)^2 .* Tstar[:].^2 .* (1 .- grounded_fraction[:])
        
    elseif ~(qmr.melt_partial_cell) && (qmr.flocal) #no partial cell melting and local
        basal_melt[grounded_fraction .== 0] .=  qmr.γT .* (qmr.ρw * qmr.c / qmr.ρi /qmr.L)^2 .* Tstar[grounded_fraction .== 0].^2
        basal_melt[.~(grounded_fraction .== 0)] .= 0 

    elseif (qmr.melt_partial_cell) && ~(qmr.flocal) #partial cell melting and non-local 
        basal_melt[:] .=  qmr.γT .* (qmr.ρw * qmr.c / qmr.ρi /qmr.L)^2 .* Tstar[:].^2 .* (1 .- grounded_fraction[:])

    elseif ~(qmr.melt_partial_cell) && ~(qmr.flocal) #no partial cell and non-local
        basal_melt[grounded_fraction .== 0] .=  qmr.γT .* (qmr.ρw * qmr.c / qmr.ρi /qmr.L)^2 .* Tstar[grounded_fraction .== 0] .* Tstar_shelf_mean
        basal_melt[.~(grounded_fraction .== 0)] .= 0 
    end
    basal_melt[:] .= basal_melt[:].* 365.25*24*60*60
    return nothing
end
