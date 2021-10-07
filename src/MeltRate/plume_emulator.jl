struct PlumeEmulator{T <: Real} <: AbstractMeltRate
    α :: T         #calibration coefficient 
    λ1 :: T        #liquidus slope 
    λ2 :: T        #liquidus intercept 
    λ3 :: T        #liquidus pressure coefficient
    E0 :: T        #entrainment coefficient
    Cd :: T        #drag coeffecient
    Γ_TS :: T      #combined Stanton number
    L :: T         #latent heat of fusion for ice 
    c :: T         #specific heat capacity of ocean water 
    βs :: T        #haline contraction coefficient 
    βt :: T        #thermal expansion coefficient 
    g :: T         #gravity
    ρi :: T        #ice density 
    ρw :: T        #water density
    Ta :: Function #Function specifying ambient temperature  
    Sa :: Function #Function specifying ambient salinity 
end

"""
    PlumeEmulator(; <kwargs>)

Construct a plume emulator to prescribe the melt rate in WAVI

Keyword arguments
=================
- α: calibration coefficient: constant prefactor used to calibrate melt rate. Defaults to 0.73 according to Favier 2019 (10.5194/gmd-12-2255-2019)
- λ1: liquidus slope (units ∘C)
- λ2: liquidus intercept (units ∘C)
- λ3: liquidus pressure coefficient (units K/m)
- E0: entrainment coefficient (dimensionless)
- Cd: drag coefficient (dimensionless)
- Γ_TS: combined Stanton number (dimensionless)
- L: latent heat of fusion for ice (units: J/kg)
- c: specific heat capacity of water (units J/kg/K)
- βs: haline contraction coefficient (units psu^(-1))
- βt: thermal expansion coefficient (units C^{-1})
- g: gravitational acceleration (units m / s^2)
- ρi: ice density (units kg/m^3)
- ρw: water density (units kg/m^3)
- Ta: ambient temperature function, defaults to the ISOMIP warm0 scenario. Ta must be a function of a single variable (depth) [time dependence not yet implemented in WAVI.jl]
- Sa: ambient salnity function, defaults to the ISOMIP warm0 scenario.  Sa must be a function of a single variable (depth) [time dependence not yet implemented in WAVI.jl]
"""
function PlumeEmulator(;
                        α = 0.73,
                        λ1 = -5.73e-2,
                        λ2 = 8.32e-4,
                        λ3 = 7.61e-4,
                        E0 = 3.6e-2,
                        Cd = 2.5e-3,
                        Γ_TS = 0.0118,
                        L = 3.35e5,
                        c = 3.974e3, 
                        βs = 7.86e-4,
                        βt = 3.87e-5,
                        g = 9.81,
                        ρi = 918.0,
                        ρw = 1028.0,
                        Ta = isomip_warm0_temp,
                        Sa = isomip_warm0_salinity)

    #check that Ta and Sa accept a single argument
    try Ta(0.0)
    catch
        ArgumentError("Ambient temperature function Ta must accept a single argument")
    end
    try Sa(0.0)
    catch
        ArgumentError("Ambient salinity function Sa must accept a single argument")
    end

    return PlumeEmulator(α, λ1, λ2, λ3,E0, Cd, Γ_TS, L, c, βs, βt,g, ρi, ρw, Ta, Sa)
    
end


function update_melt_rate!(plume_emulator::PlumeEmulator, fields, grid)
    @unpack basal_melt, h, b, grounded_fraction = fields.gh #get the ice thickness and grounded fraction
    set_plume_emulator_melt_rate!(basal_melt, 
                                h, 
                                grounded_fraction, 
                                b,
                                grid.dx,
                                grid.dy,
                                plume_emulator)
    return nothing
end


"""
    set_plume_emulator_melt_rate(melt, h, grounded_fraction)

Set the melt rate in according to Lazeroms 2018. All plume parameters passed in the plume model emulator (pme)
"""
function set_plume_emulator_melt_rate!(melt, h, grounded_fraction,bathy, dx, dy,pme)
    nx, ny = size(h)
    zbf = @. -(pme.ρi/pme.ρw)*h.*(1-grounded_fraction) .+ grounded_fraction*bathy #ice draft if floating everywhere
    ∂zb∂x, ∂zb∂y = get_slope(zbf,dx, dy) #returns the partial derivates of base in both directions

    #loop over each grid point
    for i = 1:nx
        for j = 1:ny
            melt[i,j] =  pme_pointwise_melt(grounded_fraction,zbf, ∂zb∂x,∂zb∂y, pme, i, j)
        end 
    end 

    return melt
end

"""
    pme_pointwise_melt(grounded_fraction,zbf, ∂zb∂x,∂zb∂y, pme, i, j)

Return the pme prediction of melt rate at single grid point (i,j). 
____________
Arguments:
- grounded_fraction: matrix of grounded fraction at grid points 
- zbf: ice shelf basal elevation
- ∂zb∂x: slope of ice shelf draft at grid pts in x
- ∂zb∂y: slope of ice shelf draft at grid pts in y
- pme: plume model emulator structure. Holds all parameters. 
- i,j: specific grid pt
"""
function pme_pointwise_melt(grounded_fraction,zbf, ∂zb∂x,∂zb∂y, pme, i, j)
    if grounded_fraction[i,j] == 0. #i.e. if its floating
        #initialize directions 
        dirs = get_plume_directions() #array of directions
        nd = size(dirs)[1] #number of search directions
        norm_dirs = dirs ./ norm.(eachrow(dirs)) #normalized directions    

        #initialize slope and gl depth
        zgl = zeros(1, nd) #grounding line depth
        s = zeros(1, nd) #basal slopes
        iscount = trues(1, nd) #do these points count or not?
        iscount_s = trues(1,nd)
        iscount_zgl = trues(1,nd)

        #loop over directions
        for k = 1:size(dirs)[1]
            #compute the slope in this direction
            s[k] = -dot(norm_dirs[k,:], [∂zb∂x[i,j] ∂zb∂y[i,j]]) #note negative sign
            if s[k] <= 0
                iscount[k] = false
                iscount_s[k] = false
            end #[this syntax count be tightened]

            #compute gl depth in this direction
            zgl[k] = get_pme_gl_depth(dirs[k,:], grounded_fraction, zbf, i, j)
            if isnan(zgl[k]) #no gl point found in this direction
                iscount[k] = false
                iscount_zgl[k] = false
            elseif zgl[k] >= zbf[i,j] #if gronding line higher than shelf point
                iscount[k] = false 
                iscount_zgl[k] = false
            end
        end

        #compute effective gl depth and angle
        N = sum(iscount)
        zgl_eff = sum(zgl[iscount])/N;
        α = atan(sum(s[iscount])/N);

        #update the melt rate
        if N > 0
            Sa_gl = pme.Sa(zgl_eff)
            Sa_local = pme.Sa(zbf[i,j])
            Ta_local = pme.Ta(zbf[i,j])
            Tf0 = pme.λ1*Sa_gl + pme.λ2 + pme.λ3 * zgl_eff
            M0 = lazeroms_melt_rate_prefactor(pme.βs, pme.βt, pme.g, pme.Cd, pme.Γ_TS, pme.L, pme.c, pme.λ1, pme.λ3, pme.E0, α, Sa_local, Ta_local)
            x = lazeroms_dimensionless_distance(pme.λ3, zbf[i,j], zgl_eff, Ta_local, Tf0)
            melt = M0 * pme.α * (Ta_local - Tf0)^2 * lazeroms_dimensionless_melt(x) * 365.25* 24 * 60^2
        else
            melt = 0.
        end
    else
        melt = 0.
    end
    return melt #, zgl, s, iscount, iscount_s, iscount_zgl
end


"""
    get_slope(h, dx, dy)

Returns the partial derivatives of h with respect to x and y as matrices. Uses centred FD at interior points and one sided FD at edge points (both second order)

"""
function get_slope(h, dx, dy)
    ∂h∂x = zeros(size(h))
    ∂h∂x[2:end-1,:] = (-h[1:end-2,:] + h[3:end,:])/ 2/dx 
    ∂h∂x[1,:] = (-3/2 *h[1,:] + 2*h[2,:] - 1/2 * h[3,:] )/dx
    ∂h∂x[end,:] = (1/2 * h[end-2, :] - 2 *h[end-1,:] + 3/2 * h[end,:])/dx

    ∂h∂y = zeros(size(h))
    ∂h∂y[:,2:end-1]= (-h[:,1:end-2] + h[:,3:end])/2 /dy
    ∂h∂y[:,1] = (-3/2 * h[:,1] + 2 * h[:,2] - 1/2 * h[:,3])/dy
    ∂h∂y[:,end] = (1/2 * h[:,end-2] - 2*h[:,end-1] + 3/2 * h[:,end] )/dy;

    return ∂h∂x, ∂h∂y
end

function get_plume_directions()
    directions = [0 1; #N
                  1 2; #NNE 
                  1 1; #NE 
                  2 1; #ENE 
                  1 0; #E 
                  2 -1; #ESE 
                  1 -1; #SE 
                  1 -2; #SSE 
                  0 -1; #S 
                  -1 -2; #SSW 
                  -1 -1; #SW 
                  -2 -1; #WSW 
                  -1 0; #W 
                  -2 1; #WNW 
                  -1 1; #NW 
                  -1 2] #NNW 
    return directions
end

function lazeroms_melt_rate_prefactor(βs,βt, g, Cd, Γ_TS, L, c,λ1, λ3,  E0, α, Sa, Ta)
    cp1 = (L/c)/sqrt(Cd)/Γ_TS * βt/βs / Sa
    cp2 = -λ1 * βt/βs
    ct = cp2/cp1
    M0 = sqrt(βs * Sa * g) / sqrt(λ3 * (L/c)^3) *
         ((1 -  cp1 * sqrt(Cd)*Γ_TS)/(Cd + E0 * sin(α)))^(1/2) *
         (sqrt(Cd) * Γ_TS * E0 * sin(α)/ (sqrt(Cd) * Γ_TS + ct + E0 * sin(α)))^(3/2)
    return M0
end

"""
    lazeroms_dimensionless_distance(λ3, zb, zgl, Ta, Tf0)

The dimensionless distance in the Lazeroms plume emulator formulation (see equation 28b in Lazeroms et al. 2019)
Note that we include only the first term in the expansion (i.e. no contribution from term C_ϵ)
"""
lazeroms_dimensionless_distance(λ3, zb, zgl, Ta, Tf0) = λ3 * (zb - zgl) / (Ta - Tf0)

"""
    lazeroms_dimensionless_melt(x)

Dimensionless melt function, equation 26 in Lazeroms et al. 2019 
"""
lazeroms_dimensionless_melt(x) = (3*(1-x)^(4/3) - 1)*(1- (1-x)^(4/3))^(1/2) / 2 / sqrt(2)


"""
    get_pme_gl_depth(direction, grounded_fraction, zbf, i,j)

Get the effective grounding line depth at shelf point i,j in the specified direction.
Note that this function does not yet apply the interpolation described in Lazeroms et al. 2018
"""
function get_pme_gl_depth(direction, grounded_fraction, zbf, i,j)
    nx, ny = size(grounded_fraction)

    #initialize
    zgl = NaN
    isterm_gl = false #have we found a grounded point?
    isterm_domain = false #have we gone outsdie of the domain?
    p = i 
    q = j
    
    #workflow: (1) check we aren't grounded (before step as may be grounded to begin with)
    #          (2) make a step
    #          (3) check we haven't gone out of the domain
    #          (4) loop if still in domain and not grounded
    while ~(isterm_domain) && ~(isterm_gl)
        #check whether search point is grounded, trigger terminated if so 
        grounded_fraction[p,q] > 1 - 1e-3 ?  isterm_gl = true : isterm_gl = false
        if isterm_gl #we've found a grounded point
            zgl = zbf[p,q]
        end

        #update search point
        p = p + direction[1]
        q = q + direction[2]

        #have we gone out of the domain
        isin_domain(p,q,nx,ny) ? isterm_domain = false : isterm_domain = true #if the first step is out of the domain, terminate
    end #end while loop
    return zgl 
end


isin_domain(p,q,nx,ny) = ((p <= nx) && (p>0) && (q >0) && (q <= ny)) ? true : false


