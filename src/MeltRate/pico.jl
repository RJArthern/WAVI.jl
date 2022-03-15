struct PICO{T <: Real, N <: Int} <: AbstractMeltRate
    γT :: T                       #(calibrated) heat exchange velocity
    λ1 :: T                       #liquidus slope 
    λ2 :: T                       #liquidus intercept 
    λ3 :: T                       #liquidus pressure coefficient
    L :: T                        #latent heat of fusion for ice 
    c :: T                        #specific heat capacity of ocean water 
    βs :: T                       #haline contraction coefficient 
    βt :: T                       #thermal expansion coefficient
    g  :: T                       #gravity 
    ρi :: T                       #ice density 
    ρw :: T                       #water density
    C  :: T                       #overturning strength
    nbox :: N                     #number of boxes
    S0 :: T                       #input salinity
    T0 :: T                       #input temperature
    ice_front_mask::Array{Bool,2} #specify ice front location
    use_box_mean_depth::Bool      #flag to specify whether to use the mean depth in a box for the calculation of pk 
end

"""
    function PICO(; <kwargs>)

Construct a pico structure to prescribe the melt rate in WAVI

Keyword arguments
=================
- γT: calibrated heat exchange velocity (units m/s). Defaults to 1.44*1e-5 according to Favier 2019 for 10 boxes (10.5194/gmd-12-2255-2019)
- λ1: liquidus slope (units ∘C)
- λ2: liquidus intercept (units ∘C)
- λ3: liquidus pressure coefficient (units K/m)
- L: latent heat of fusion for ice (units: J/kg)
- c: specific heat capacity of water (units J/kg/K)
- βs: haline contraction coefficient (units psu^(-1))
- βt: thermal expansion coefficient (units ∘C^{-1})
- g: gravitational acceleration (units m / s^2)
- ρi: ice density (units kg/m^3)
- ρw: water density (units kg/m^3)
- C: overturning circulation parameter (units m^6 /s /kg)
- S0: input salinity in box 0  (units PSU)
- T0: input temperature in box 0 (units ∘C)
- ice_front_mask: mask specifying the location of the ice front
- use_box_mean_depth: flag to specify whether or not to use the mean depth of the draft in each box in the computation of pressure in box
"""
function PICO(;
                γT = 1.44*1e-5,
                λ1 = -5.73e-2,
                λ2 = 8.32e-4,
                λ3 = 7.61e-4,
                L = 3.35e5,
                c = 3.974e3, 
                βs = 7.7e-4,
                βt = 7.5e-5,
                g = 9.81, 
                ρi = 918.8,
                ρw = 1028.0,
                C = 1.0e6,
                nbox = 10,
                S0 = 34.6, 
                T0 = 0.0,
                ice_front_mask = nothing,
                use_box_mean_depth = false)

    #check that we've input an ice front
    ~(ice_front_mask === nothing) || throw(ArgumentError("Must input an ice front matrix"))

    #map the ice front mask to Boolean 
    try
        ice_front_mask = convert(Array{Bool,2}, ice_front_mask)
    catch
        throw(ArgumentError("Ice front mask must be able to be cast as a Boolean matrix"))
    end
    return PICO(γT, λ1, λ2, λ3, L, c, βs, βt, g, ρi, ρw, C, nbox, S0, T0, ice_front_mask, use_box_mean_depth)
end

"""
    function update_melt_rate!(pico::PICO, fields, grid, clock)

Update the ice model melt rate for a PICO melt rate parametrization
"""
function update_melt_rate!(pico::PICO, fields, grid, clock)
    @unpack basal_melt, h, b, grounded_fraction = fields.gh #get the ice thickness and grounded fraction
    
    #check that the ice front mask has the correct dimensions 
    (size(pico.ice_front_mask)) == (grid.nx, grid.ny) || throw(ArgumentError("Ice front mask ($(size(ice_front_mask))) size must be same as model grid ($nx x $ny)"))
   
    #compute the ice draft
    zb = fields.gh.b .* (grounded_fraction .== 1) + - pico.ρi / pico.ρw .* fields.gh.h .* (grounded_fraction .< 1)

    set_pico_melt_rate!(basal_melt, 
                        grounded_fraction,
                        pico, 
                        grid,
                        zb)
    return nothing
end

"""
    function set_pico_melt_rate!(melt, grounded_fraction, pico, grid, zb)

Update the melt rate (stored in melt) using the specified pico model.
"""
function set_pico_melt_rate!(melt,grounded_fraction, pico, grid, zb)
    #determine which box certain parts are in (i.e. are you floating, and what is the distance to the ice front?)
    grounding_line_mask    = get_gl_mask(grounded_fraction) #boolean matrix specifying where the grounding line is
    boxes = get_pico_boxes(grounding_line_mask, pico.ice_front_mask, grounded_fraction, pico.nbox)

    #set the melt rate in boxes sequentially
    T_out, S_out,q = set_box_one_melt_rate!(melt, boxes,pico, grid.dx,grid.dy, zb)
    for k = 2:pico.nbox
        T_out, S_out = set_box_k_melt_rate!(melt,k, T_out, S_out,q, boxes,pico, grid.dx,grid.dy, zb)
    end
end

"""
    function get_pico_boxes()

Return an array of boxes for the PICO melt rate. Returns 0 for grounded points.
"""
function get_pico_boxes(grounding_line_mask, ice_front_mask, grounded_fraction, nbox)
    nx, ny = size(grounding_line_mask)
    boxes = zeros(nx,ny)
    r = zeros(nx,ny)

    inds_gl = Tuple.(findall(!iszero, grounding_line_mask))
    inds_if = Tuple.(findall(!iszero, ice_front_mask))     #indices of ice front and grounding line
    for i = 1:nx
        for j = 1:ny
            if grounded_fraction[i,j] == 0 & ~(grounding_line_mask[i,j]) #if we're in the shelf and not at a grounding line point
                #compute distance to the gl
                dgl = findmin(((i .- first.(inds_gl)).^2 + (j .- last.(inds_gl)).^2).^(1/2))[1]
                #compute distance to the ice front
                dif = findmin(((i .- first.(inds_if)).^2 + (j .- last.(inds_if)).^2).^(1/2))[1]
                #assign box
                r[i,j] = dgl / (dgl + dif) #equation 10 in Reese et al 2018
                boxes[i,j] = ceil(r[i,j] * nbox)
            end
         end
    end
    return boxes
end

""" 
    function get_gl_mask(grounded_fraction)

Return a matrix with ones corresponding to locations of grounding line points
Grounding line points determines as those points with grounded fraction < 1, and adjacted to a point with grounded fraction = 1
"""
function get_gl_mask(grounded_fraction)
    nx, ny = size(grounded_fraction)
    gl_mask = falses(nx,ny) 
    for i = 1:nx
        for j = 1:ny
            if (grounded_fraction[i,j] < 1)
                #compute the box of adjacent points
                adj_box_grounded = get_adj_box(i,j,grounded_fraction)
                #if any adjacent point has grounded_fraction == 1, tag this point as being in gl
                if any(adj_box_grounded .== 1)
                    gl_mask[i,j] = true
                end
            end
        end
    end
    return gl_mask
end

"""
    function get_adj_box(i,j,m)

Returns the entries of m which are adjacent to grid point indexed by i, j, taking into account edges of matrix
"""
function get_adj_box(i,j, m)
    nx, ny = size(m)
    lxn = min(1,i-1) #length of box negative x
    lxp = min(1,nx - i) #length of box positive x
    lyn = min(1,j-1)
    lyp = min(1,ny - j)
    box = m[((i-lxn):(i+lxp)), ((j - lyn):(j+lyp))]
    return box
end

"""
    function set_box_one_melt_rate!(melt, boxes, pico, dx, dy, zb)

Set the melt rate in pico box one (adjacent to the grounding line). Box location specified by boxes array.
Returns the mean temperature T1bar and salinity S1bar in box one (used as boundary conditions on box two), and the overturning flux q (used in all subsequent boxes)
"""
function set_box_one_melt_rate!(melt,boxes,pico, dx, dy, zb)
    m = zeros(size(zb))

    #scalar quantities
    A1 = sum(sum(boxes .== 1))*dx *dy 
    g1 = A1 * pico.γT
    ν = pico.ρi / pico.ρw
    λ = pico.L / pico.c
    s  = pico.S0 / ν / λ
    if ~(pico.use_box_mean_depth) #don't use mean depth, below are vector quantities
        Tstar = pico.λ1 .* pico.S0 .+ pico.λ2 .+ pico.λ3 .* zb .- pico.T0
        c =  g1 .* Tstar / pico.C / pico.ρw / (pico.βs * s - pico.βt) #ρw = ρ* in Reese et al 2018
        b =  g1  ./2 ./ pico.C / pico.ρw ./ (pico.βs * s - pico.βt) 
        x = -b .+ sqrt.(b.^2 .- c) 
        T1 = pico.T0 .- x 
        y = pico.S0 .* x ./ ν ./ λ
        S1 = pico.S0 .- y
        m  = -pico.γT ./ ν ./ λ .* (pico.λ1 .* S1 .+ pico.λ2 .+ pico.λ3 .* zb .- T1)*365.25*24*60^2 #melt rate if all domain in box one

        #set melt rate at box one points
        melt[boxes .== 1] .= m[boxes .== 1]

        #compute the mean temperature and salinity
        T1bar = sum(T1[boxes .== 1])/length(T1[boxes .== 1])
        S1bar = sum(S1[boxes .== 1])/length(S1[boxes .== 1])

    else #use the mean depth, then below are scalar quantities
        zbf = zb
        zb = sum(zbf[boxes .== 1])/length(zbf[boxes .== 1]) #use mean depth in the box
        Tstar = pico.λ1 .* pico.S0 .+ pico.λ2 .+ pico.λ3 .* zb .- pico.T0
        c =  g1 .* Tstar / pico.C / pico.ρw / (pico.βs * s - pico.βt) #ρw = ρ* in Reese et al 2018
        b =  g1  ./2 ./ pico.C / pico.ρw ./ (pico.βs * s - pico.βt) 
        x = -b .+ sqrt.(b.^2 .- c) 
        T1 = pico.T0 .- x 
        y = pico.S0 .* x ./ ν ./ λ
        S1 = pico.S0 .- y

        #set melt rate at box one points
        m  = -pico.γT ./ ν ./ λ .* (pico.λ1 .* S1 .+ pico.λ2 .+ pico.λ3 .* zb .- T1)*365.25*24*60^2 #melt rate if all domain in box one
        melt[boxes .== 1] .= m

        #compute the mean
        T1bar = T1
        S1bar = S1
    end

    #compute the overturning
    q = pico.C * pico.ρw * (pico.βs *(pico.S0 - S1bar) - pico.βt *(pico.T0 - T1bar))
    
    return T1bar, S1bar, q
end


"""
    function set_box_k_melt_rate!(melt,k, T_prev, S_prev,q, boxes, pico, dx, dy, zb)

Set the melt rate in pico box k > 1.  (Box location specified by boxes array.)
Returns the mean temperature and salinity in this box (used as boundary conditions on the next box)
"""
function set_box_k_melt_rate!(melt,k, T_prev, S_prev,q,boxes,pico, dx, dy, zb)
    m = zeros(size(zb))

    #scalar quantities
    Ak = sum(sum(boxes .== k))*dx *dy 
    g1 = Ak * pico.γT
    ν = pico.ρi / pico.ρw
    λ = pico.L / pico.c
  
    #compute salinity and temp in this box
    if ~(pico.use_box_mean_depth) #don't use mean depth, below are vector quantities
        Tstar = pico.λ1 .* S_prev .+ pico.λ2 .+ pico.λ3 .* zb .- T_prev
        x = -g1 .* Tstar ./ (q + g1 - g1*pico.λ1 * S_prev / ν / λ)
        y = S_prev .* x ./ ν / λ
        Tk = T_prev .- x
        Sk = S_prev .- y
    
        # compute and set melt rate in this box
        m  = -pico.γT ./ ν ./ λ .* (pico.λ1 .* Sk .+ pico.λ2 .+ pico.λ3 .* zb .- Tk)*365.25*24*60^2 #melt rate if all domain in box one
        melt[boxes .== k] .= m[boxes .== k]

        #update the mean temperature and salinity for the next box
        T_out = sum(Tk[boxes .== k])/length(Tk[boxes .== k])
        S_out = sum(Sk[boxes .== k])/length(Sk[boxes .== k])

    else #use the mean depth, then below are scalar quantities
        zbf = zb
        zb = sum(zbf[boxes .== k])/length(zbf[boxes .== k]) #use mean depth in the box
        Tstar = pico.λ1 .* S_prev .+ pico.λ2 .+ pico.λ3 .* zb .- T_prev
        x = -g1 .* Tstar ./ (q + g1 - g1*pico.λ1 * S_prev / ν / λ)
        y = S_prev .* x ./ ν / λ
        Tk = T_prev .- x
        Sk = S_prev .- y
    
        # compute and set melt rate in this box
        m  = -pico.γT ./ ν ./ λ .* (pico.λ1 .* Sk .+ pico.λ2 .+ pico.λ3 .* zb .- Tk)*365.25*24*60^2 #melt rate if all domain in box one
        melt[boxes .== k] .= m

        T_out = Tk
        S_out = Sk
    end    
    return T_out, S_out
end