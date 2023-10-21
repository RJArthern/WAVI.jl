using Statistics

struct IdealizedAnthroMeltRate{T <: Real, N <: Integer} <: AbstractMeltRate
    bump_width :: T               #lengthscale of the quadratic bump in forcing
    bump_time  :: T               #time of the center of the bump
    bump_amplitude :: T           #amplitude of the bump
    per_century_trend :: T        #per century trend in the forcing
    trend_onset :: T              #time of onset of the trend
    pc_max :: T                   #maximum pycnocline position (without a trend) 
    pc_min :: T                   #minimum pycncoline position (without a trend)
    M :: T                        #melt rate prefactor = γT .* (ρw * c / ρi /L)^2
    λ1 :: T                       #liquidus slope 
    λ2 :: T                       #liquidus intercept 
    λ3 :: T                       #liquidus pressure coefficient
    melt_partial_cell :: Bool     #Flag for melting applied to partial cells or not
    random_seed :: N              #(integer) random seed for the stochastic process
    rf_threshold :: T             #threshold for the autoregressive process to be capped at
    r :: T                        #autocorrelation of AR1 process
    Sl :: T                       #salinity of the lower layer 
    Su :: T                       #salinity of the upper layer
    Tl :: T                       #temperature of the upper layer
    Tu :: T                       #temperature of the lower layer
    pw :: T                       #width of the pycnocline
end


"""
    IdealizedAnthroMeltRate(; ,kwargs)

Construct a melt rate object for the idealized anthropogenic forcing experiments. This is a quadratic dependency of melting on depth with a variable prefactor, trend in forcing, and size of bump

Keyword arguments 
=================
- bump_width :: T               #lengthscale of the quadratic bump in forcing
- bump_time  :: T               #time of the center of the bump
- bump_amplitude :: T           #amplitude of the bump
- per_century_trend :: T        #per century trend in the forcing
- trend_onset :: T              #time of onset of the trend
- pc_max :: T                   #maximum pycnocline position (without a trend) 
- pc_min :: T                   #minimum pycncoline position (without a trend)
- M :: T                        #melt rate prefactor = γT .* (ρw * c / ρi /L)^2
- λ1 :: T                       #liquidus slope 
- λ2 :: T                       #liquidus intercept 
- λ3 :: T                       #liquidus pressure coefficient
- melt_partial_cell             #Flag for melting applied to partial cells or not
- random_seed                   #random seed for the stochastic process
- rf_threshold                  #threshold for the autoregressive process to be capped at
- r                             #autocorrelation of AR1 process
- Sl                            #salinity of the lower layer 
- Su                            #salinity of the upper layer
- Tl                            #temperature of the upper layer
- Tu                            #temperature of the lower layer
- pw                            #width of the pycnocline
"""

function IdealizedAnthroMeltRate(;
                                bump_width = 1.0e-5,
                                bump_time  = 0.0,
                                bump_amplitude = 0.0, 
                                per_century_trend = 0.0, 
                                trend_onset = 0.0, 
                                pc_max = -400.0,
                                pc_min = -600.0,
                                M = 0.56, #value from Bradley et al. 2023
                                λ1 = -5.73e-2,
                                λ2 = 8.32e-4,
                                λ3 = 7.61e-4,
                                melt_partial_cell = true, 
                                random_seed = 1, 
                                rf_threshold = 1.5, 
                                r = 0.8, 
                                Sl = 34.6,
                                Su = 34.0,
                                Tl = 1.0, 
                                Tu = -1.2,
                                pw = 400.0)

    return IdealizedAnthroMeltRate(bump_width, bump_time, bump_amplitude, per_century_trend, trend_onset,
                    pc_max, pc_min, M, λ1, λ2,λ3, melt_partial_cell, random_seed,rf_threshold, r ,
                    Sl, Su, Tl, Tu, pw)
end


"""
    function update_melt_rate!(idealized_anthro_melt_rate, fields, grid, clock)

Wrapper script to update the basal melt rate for an idealized antho melt rate object.

"""
function update_melt_rate!(idealized_anthro_melt_rate::IdealizedAnthroMeltRate, fields, grid, clock)
    @unpack basal_melt, h, b, grounded_fraction = fields.gh #get the ice thickness and grounded fraction
 
    #compute the ice draft
    ρi = 918.0
    ρw = 1028.0
    zb = fields.gh.b .* (grounded_fraction .== 1) + - ρi / ρw .* fields.gh.h .* (grounded_fraction .< 1)

    set_idealized_anthro_melt_rate!(basal_melt,
                           idealized_anthro_melt_rate,
                            zb, 
                            grounded_fraction, clock.time)
    return nothing
end

"""
function set_idealized_anthro_melt_rate!(basal_melt, idealized_anthro_melt_rate, zb,  grounded_fraction, t)
    
    Set the melt rate field (basal_melt) for the idealized_anthro_melt_rate object at time t

"""
function set_idealized_anthro_melt_rate!(basal_melt,
                           idealized_anthro_melt_rate,
                            zb, 
                            grounded_fraction, t)       

    #get the current pycncoline position
    pc = pc_position(t,
                idealized_anthro_melt_rate.r, 
                idealized_anthro_melt_rate.random_seed, 
                idealized_anthro_melt_rate.rf_threshold, 
                idealized_anthro_melt_rate.pc_max, 
                idealized_anthro_melt_rate.pc_min, 
                idealized_anthro_melt_rate.per_century_trend, 
                idealized_anthro_melt_rate.trend_onset,
                idealized_anthro_melt_rate.bump_width, 
                idealized_anthro_melt_rate.bump_amplitude, 
                idealized_anthro_melt_rate.bump_time) 
        
    #compute the salinity and temperature
    Sa_shelf = get_Sa.(zb,
                idealized_anthro_melt_rate.Sl,
                idealized_anthro_melt_rate.Su,
                pc,
                idealized_anthro_melt_rate.pw) #basal ambient salinity at all points (restrict melting to only floating points below)

    Ta_shelf = get_Ta.(zb,
                idealized_anthro_melt_rate.Tl,
                idealized_anthro_melt_rate.Tu,
                pc,
                idealized_anthro_melt_rate.pw) #basal ambient temperature at all points
    
    Tstar = idealized_anthro_melt_rate.λ1 .* Sa_shelf .+ idealized_anthro_melt_rate.λ2 .+ idealized_anthro_melt_rate.λ3 .* zb .- Ta_shelf #local freezing temperature 

    # set the melt rate
    if (idealized_anthro_melt_rate.melt_partial_cell) #we do have partial cell melting 
        basal_melt[:] .= idealized_anthro_melt_rate.M .* Tstar[:].^2 .* (1 .- grounded_fraction[:])
    else  #no partial cell melting an
        basal_melt[grounded_fraction .== 0] .=  idealized_anthro_melt_rate.M .* Tstar[grounded_fraction .== 0].^2 #set melt rate on fully floating points 
        basal_melt[.~(grounded_fraction .== 0)] .= 0 #check that anywhere with grounded_fraction not equal to 1 gets no melting

    end

    return nothing
end

"""
    function generate_random_forcing_anomaly(t, r, rs)

Generate values of the random forcing as an AR1 process with autocorrelation r at time point t

"""
function generate_random_forcing_anomaly(t, r, random_seed)
    N = 1000  # number of timeseries in the series to be sampled(can be big)
    dt = 1    # script assumes this is 1, no guarantees it is robust if changed!
    td = 0:dt:(N-1)*dt

    # frequency domain
    df = 1 / (dt * N)
    f0 = 0.5 / dt
    f1 = 0:df:f0

    # persistence parameters
    beta = 0.5
    tau = dt / (1 - r)
    p0 = 1

    # power spectra zz(analytical)
    pb = sqrt.(p0 .* (f0 ./ f1[2:end]) .^ beta)
    pr = sqrt.(p0 ./ (1 .+ r^2 .- 2r .* cos.(2π .* dt .* f1[2:end])))

    # random phases - must be conjugate symmetric
    phase = zeros(Complex{Float64}, N)
    rng = MersenneTwister(random_seed)
   # rng = seed!(MersenneTwister(random_seed))
    pos_freq = pb[2:div(N, 2)]
    ranphase = Random.rand(rng, size(pos_freq)[1])
    phase[2:length(pos_freq)+1] .= 1im .* 2π .* ranphase
    phase[N-length(pos_freq)+1:end] .= conj(reverse(phase[2:length(pos_freq)+1]))

    pb2 = [pb; reverse(pb[1:div(N, 2)])]
    pr2 = [pr; reverse(pr[1:div(N, 2)])]

    # combine phase and power spectra
    pbrand = pb2 .* exp.(phase)
    prrand = pr2 .* exp.(phase)

    ar1 = real(AbstractFFTs.ifft(prrand))
    ar1 = ar1.-mean(ar1);
    ar1 = ar1./Statistics.std(ar1);

    # interpolate to target value
    linterp = Interpolations.LinearInterpolation(td, ar1)
    rf = linterp(t)
 
    return rf
end

"""
    function get_random_pc_component()

Map the value of the autoregressive process to the pycnocline position, i.e. get the internal part of the pycnoline position
"""
function get_random_pc_component(t, r, random_seed, rf_threshold, pc_max, pc_min)
    rf = generate_random_forcing_anomaly(t, r, random_seed)
    if rf  > rf_threshold
        rf = rf_threshold
    end

    if rf < -rf_threshold
        rf = -rf_threshold
    end
    pc_pos = (pc_min+pc_max)/2 +  (pc_max - pc_min)*rf/2/abs(rf_threshold)

    return pc_pos
end

"""
    function get_trend_pc_component(t, per_century_trend, trend_onset)

Return the trend component of the pycncoline position at a given t.
"""
function get_trend_pc_component(t, per_century_trend, trend_onset)
    pc_pos = 0
    if t < trend_onset
        pc_pos = 0
    else
        pc_pos = per_century_trend/100 * (t - trend_onset)
    end

    return pc_pos 
end

get_bump_pc_component(t, bump_width, bump_amplitude, bump_time) = bump_amplitude*exp(-(t - bump_time)^2 /2 /bump_width^2)



pc_position(t,r, random_seed, rf_threshold, pc_max, pc_min, per_century_trend, trend_onset,bump_width, bump_amplitude, bump_time) = get_random_pc_component(t,r, random_seed, rf_threshold, pc_max, pc_min) + get_trend_pc_component(t, per_century_trend, trend_onset) + get_bump_pc_component(t, bump_width, bump_amplitude, bump_time)

"""
    function get_Ta(z,Tl,Tu,pc)

return the ambient temperature as a function of z for a given Tl (lower layer temperature), Tu (upper layer temperature), and pycncoline centre position pc and pycncoline width pw

"""
function get_Ta(z,Tl,Tu,pc, pw)
    Ta = Tu
    if z > (pc + pw)
        Ta = Tu

    elseif z < (pc - pw)
        Ta = Tl

    else
        Ta = Tl +(z - (pc-pw))*(Tu - Tl)/(2*pw)

    end

    return Ta
end

"""
    function get_Sa(z,Tl,Tu,pc)

return the ambient salinity as a function of z for a given Sl (lower layer temperature), Su (upper layer temperature), and pycncoline centre position pc and pycncoline width pw

"""
function get_Sa(z,Sl,Su,pc, pw)
    Sa = Su
    if z > (pc + pw)
        Sa = Su

    elseif z < (pc - pw)
        Sa = Sl

    else
        Sa = Sl +(z - (pc-pw))*(Su - Sl)/(2*pw)

    end

    return Sa
end