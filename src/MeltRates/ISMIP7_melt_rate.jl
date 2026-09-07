export ISMIP7MeltRate

using WAVI: AbstractClimateForcing
using WAVI.ClimateForcing
using NCDatasets 

struct ISMIP7MeltRate{P <: String, V<: String, T <: Real} <: AbstractMeltRate
    so_prefix::P                       #Pass the string which is the prefix of the salinity forcing (e.g. "so_AIS_CESM2-WACCM_ssp585_ocean_v3_16000m_")
    tf_prefix::P                       #Pass the string which is the prefix of the thermal forcing (e.g. "tf_AIS_CESM2-WACCM_ssp585_ocean_v3_16000m_")
    thetao_prefix::P                   #Pass the string which is the prefix of the thetao forcing (e.g. "thetao_AIS_CESM2-WACCM_ssp585_ocean_v3_16000m_")
    so_varname::V                      #Variable name for the salinity in the forcing file (defualt: "so")
    tf_varname::V                      #Variable name for the thermal forcing in the forcing file (defualt: "tf")
    thetao_varname::V                  #Variable name for the temperature in the forcing file (defualt: "thetao")
    K::T                                #tuning parameter (see here for info: https://drive.google.com/file/d/1SygMQte-7XgKj4e-Hpyj1tHi6XguShCP/view and https://tc.copernicus.org/articles/16/4931/2022/tc-16-4931-2022.pdf)
    shelf_slope::T                      #average slope of ice shelf base (can be local in ISMIP7, but we're doing with the average)
    ρ_ocean::T                          #density of the ocean
    ρ_ice::T                            #density of the ice 
    c_ocean::T                             #thermal heat paracity of the ocean 
    L_ice::T                            #ice latent heat of fusion 
    β_s::T                              #Salt contraction coefficient
    g::T                                #gravity 
    f::T                                #Coriolis parameter
    S_loc::Union{Array{T,3}, Nothing}   #array to hold the 3d salinity, to be read in from ISMIP7 forcing
    T_loc::Union{Array{T,3}, Nothing}   #array to hold the 3d temperature, to be read in from ISMIP7 forcing
    Tf_loc::Union{Array{T,3}, Nothing}  #array to hold the 3d thermal forcing, to be read in from ISMIP7 forcing
    z_forcing::Union{Array{T,1}, Nothing} #array to hold the z co-ordinates from the ISMIP7 forcing
    melt_partial_cell::Bool             #Flag for melting applied to partial cells or not
    path_to_forcing :: String
    # Global-grid index ranges.
    # Ocean fields stay on the full forcing grid; these map local (i,j) to global.
    x_indices::Union{Nothing, UnitRange{Int}}
    y_indices::Union{Nothing, UnitRange{Int}}
end

"""
    ISMIP7MeltRate(;, kwargs)

Construct an ISMIP7MeltRate object to prescribe the melt rate in WAVI

Keyword arguments
=================
- climate_forcing: decriptor of climate forcing (ISMIP7_ANOMALY, ISMIP7_OCX, or ISMIP7_CONTROL)
- K: melt rate tuning parameter (see here for info: https://drive.google.com/file/d/1SygMQte-7XgKj4e-Hpyj1tHi6XguShCP/view and https://tc.copernicus.org/articles/16/4931/2022/tc-16-4931-2022.pdf)
- shelf_slope: average slope of ice shelf base (can be local in ISMIP7, but we're doing with the average)
- ρ_ocean: density of the ocean water
- ρ_ice: density of the ice 
- c_ocean: thermal heat paracity of the ocean 
- L_ice: ice latent heat of fusion 
- Salt: contraction coefficient
- g: gravity 
- f: Coriolis parameter
- S_loc: array to hold the 3d salinity, to be read in from ISMIP7 forcing
- array: to hold the 3d temperature, to be read in from ISMIP7 forcing
- Tf_loc: array to hold the 3d thermal forcing, to be read in from ISMIP7 forcing
- z_forcing: array to hold the z co-ordinates from the ISMIP7 forcing
- melt_partial_cell: Flag for melting applied to partial cells or not
- path_to_forcing: path to the forcing files
- x_indices: Maps local field x indices to global NetCDF x indices when on a subdomain.
- y_indices: Maps local field y indices to global NetCDF y indices when on a subdomain.
"""
function ISMIP7MeltRate(; 
                        climate_forcing = nothing,
                        K = 1.0e-5,
                        shelf_slope = 1.e-3,
                        ρ_ocean = 1028.0,
                        ρ_ice = 918.0,
                        c_ocean = 3974.0,
                        L_ice = 3.34e5,
                        β_s = 7.7e-4,
                        g = 9.81,
                        f = -1.0e-4,
                        S_loc = nothing,
                        T_loc = nothing,
                        Tf_loc = nothing, 
                        z_forcing = nothing,
                        melt_partial_cell = false, 
                        path_to_forcing = "./",
                        x_indices = nothing,
                        y_indices = nothing)

    ~(climate_forcing === nothing) || throw(ArgumentError("You must pass a climate_forcing description"))

    if isnothing(z_forcing)
        #read in the levels
        z_levels_path = joinpath(path_to_forcing, "ocean_z_levels.nc")
        z_levels_ncfile   = NCDataset(z_levels_path)
        z_forcing = replace(z_levels_ncfile["z"][:] , missing => NaN)
    end
    
    return   ISMIP7MeltRate(climate_forcing, K, shelf_slope,ρ_ocean, ρ_ice,c_ocean, L_ice,β_s, g,f,S_loc, T_loc, Tf_loc, z_forcing, melt_partial_cell, path_to_forcing,x_indices, y_indices)
end

#Convenience Constructors
function ISMIP7MeltRate(climate_forcing::T, K, shelf_slope,ρ_ocean, ρ_ice,c_ocean, L_ice,β_s, g,f,S_loc, 
                        T_loc, Tf_loc, z_forcing, melt_partial_cell, 
                        path_to_forcing,x_indices, y_indices) where {T <: Union{ISMIP7_ANOMALY,ISMIP7_OCX,ISMIP7_CONTROL}}

    so_prefix = "so_"
    tf_prefix = "tf_"
    thetao_prefix = "thetao_"
    so_varname = "so"
    tf_varname = "tf"
    thetao_varname = "thetao"

    return ISMIP7MeltRate(so_prefix, tf_prefix, thetao_prefix, so_varname, tf_varname, thetao_varname, K, shelf_slope,ρ_ocean, ρ_ice,c_ocean, L_ice,β_s, g,f,S_loc, T_loc, Tf_loc, z_forcing, melt_partial_cell, path_to_forcing,x_indices, y_indices)
end


function update_shelf_melt_rate!(ISMIP7_melt_rate::ISMIP7MeltRate, fields, grid, clock)

    @unpack b, shelf_basal_melt, grounded_fraction, h = fields.gh 
    @unpack K, shelf_slope, ρ_ocean, ρ_ice, c_ocean, L_ice, β_s, g, f, S_loc, T_loc, Tf_loc, z_forcing, melt_partial_cell, x_indices, y_indices = ISMIP7_melt_rate

    # Ocean lookup uses host 3D arrays and scalar indexing; compute on host then copy back.
    b_h = _host(b)
    h_h = _host(h)
    gf_h = _host(grounded_fraction)

    #compute the ice draft
    zb = b_h .* (gf_h .== 1) + - ρ_ice / ρ_ocean .* h_h .* (gf_h .< 1)

    # Local draft indices into (possibly global) ocean fields
    nx, ny = size(zb)
    #compute the local thermal forcing and salinity from the climate forcing and ice draft
    idx = [argmin(abs.(z_forcing .- d)) for d in zb] #gives an nx * ny array indices which are closest depth to z in the forcing.
    ocean_is_local = size(S_loc, 1) == nx && size(S_loc, 2) == ny
    if ocean_is_local
        Tf_local_shelf = [Tf_loc[i, j, idx[i,j]] for i in 1:nx, j in 1:ny] #fills the array of thermal forcing
        S_local_shelf = [S_loc[i, j, idx[i,j]] for i in 1:nx, j in 1:ny] #fills the array of thermal forcing
    else
        is = isnothing(x_indices) ? (1:nx) : x_indices
        js = isnothing(y_indices) ? (1:ny) : y_indices
        (length(is) == nx && length(js) == ny) || throw(DimensionMismatch(
            "ISMIP7 melt x/y index ranges ($(length(is)),$(length(js))) do not match local draft size ($nx,$ny)"
        ))
        Tf_local_shelf = [Tf_loc[is[i], js[j], idx[i,j]] for i in 1:nx, j in 1:ny] #fills the array of thermal forcing
        S_local_shelf = [S_loc[is[i], js[j], idx[i,j]] for i in 1:nx, j in 1:ny] #fills the array of thermal forcing
    end


    #set the shelf melt rate
    secs_per_year = 365.25*24*60^2
    melt_h = similar(b_h)
    if melt_partial_cell
        melt_h[:] .= K * secs_per_year * shelf_slope * ρ_ocean / ρ_ice * (c_ocean/L_ice)^2 * β_s * 0.5 * g /abs(f) * S_local_shelf[:] .* abs.(Tf_local_shelf[:]) .* Tf_local_shelf[:] .*  (1 .- gf_h[:])

    elseif ~(melt_partial_cell)
        melt_h .= 0
        melt_h[gf_h .== 0] .= K * secs_per_year * shelf_slope * ρ_ocean / ρ_ice * (c_ocean/L_ice)^2 * β_s * 0.5 * g /abs(f) * S_local_shelf[gf_h .== 0] .* abs.(Tf_local_shelf[gf_h .== 0]) .* Tf_local_shelf[gf_h .== 0]
    end
    copy_onto!(shelf_basal_melt, melt_h)


    return nothing
end

function update_climate_forcing!(ISMIP7_melt_rate::ISMIP7MeltRate, grid::Grid, clock::Clock)
    @unpack S_loc,
        T_loc,
        Tf_loc,
        path_to_forcing,
        so_prefix,
        tf_prefix,
        thetao_prefix,
        so_varname,
        tf_varname,
        thetao_varname,
        x_indices,
        y_indices = ISMIP7_melt_rate

    #get the year from clock for the forcing files
    current_time = clock.time + clock.ref_time
    current_time_string = string(Int(round(current_time)))

    is, js = forcing_index_ranges(x_indices, y_indices, S_loc)

    # load in the salinity
    salinity_filename = joinpath(path_to_forcing, join([so_prefix,  current_time_string,".nc"]))
    salinity_ncfile   = NCDataset(salinity_filename)
    S_loc .= replace(salinity_ncfile[so_varname][is, js, :], missing => NaN)
    #println("read in salinity forcing file: " * salinity_filename)
    @info "read in salinity forcing file: $salinity_filename"

    # load in the temperature
    temperature_filename = joinpath(path_to_forcing, join([thetao_prefix,  current_time_string,".nc"]))
    temperature_ncfile   = NCDataset(temperature_filename)
    T_loc .= replace(temperature_ncfile[thetao_varname][is, js, :], missing => NaN)
    #println("read in temperature forcing file: " * temperature_filename)
    @info "read in temperature forcing file: $temperature_filename"


    # load in the thermal forcing
    thermal_forcing_filename = joinpath(path_to_forcing, join([tf_prefix,  current_time_string,".nc"]))
    thermal_forcing_ncfile   = NCDataset(thermal_forcing_filename)
    Tf_loc .= replace(thermal_forcing_ncfile[tf_varname][is, js, :], missing => NaN)
    #println("read in thermal-forcing forcing file: " * thermal_forcing_filename)
    @info "read in thermal-forcing forcing file: $thermal_forcing_filename 1111"


end

    
function reconstruct_on_grid(ISMIP7_melt_rate::ISMIP7MeltRate,grid::Grid)
    isnothing(ISMIP7_melt_rate.z_forcing) &&
        throw(ArgumentError("ISMIP7 melt rate needs z_forcing before reconstruct_on_grid"))
    number_of_zlevels = size(ISMIP7_melt_rate.z_forcing, 1)
    xs, ys = grid_index_ranges(ISMIP7_melt_rate.x_indices, ISMIP7_melt_rate.y_indices, grid)

    return ISMIP7MeltRate(
                ISMIP7_melt_rate.so_prefix,
                ISMIP7_melt_rate.tf_prefix,
                ISMIP7_melt_rate.thetao_prefix,
                ISMIP7_melt_rate.so_varname,
                ISMIP7_melt_rate.tf_varname, 
                ISMIP7_melt_rate.thetao_varname, 
                ISMIP7_melt_rate.K,
                ISMIP7_melt_rate.shelf_slope,
                ISMIP7_melt_rate.ρ_ocean,
                ISMIP7_melt_rate.ρ_ice,
                ISMIP7_melt_rate.c_ocean,
                ISMIP7_melt_rate.L_ice,
                ISMIP7_melt_rate.β_s,
                ISMIP7_melt_rate.g,
                ISMIP7_melt_rate.f,
                isnothing(ISMIP7_melt_rate.S_loc) ? 
                    zeros(grid.nx,grid.ny,number_of_zlevels) : 
                    ISMIP7_melt_rate.S_loc,
                isnothing(ISMIP7_melt_rate.T_loc) ? 
                    zeros(grid.nx,grid.ny,number_of_zlevels) : 
                    ISMIP7_melt_rate.T_loc,
                isnothing(ISMIP7_melt_rate.Tf_loc) ? 
                    zeros(grid.nx,grid.ny,number_of_zlevels) : 
                ISMIP7_melt_rate.Tf_loc,
                ISMIP7_melt_rate.z_forcing, 
                ISMIP7_melt_rate.melt_partial_cell,
                ISMIP7_melt_rate.path_to_forcing,
                xs,
                ys)
end

function reconstruct_on_subdomain(ISMIP7_melt_rate::ISMIP7MeltRate,grid::Grid,subdomain::NTuple{4,<: Integer}) 
    # Slice ocean arrays when they cover the ice grid, including 3D nx x ny x nz fields.
    # Keep local-to-global index maps for NetCDF hyperslab reads.
    # Compose with any existing map (e.g. Threaded tiles inside an MPI domain).
    xs, ys = subdomain_index_ranges(ISMIP7_melt_rate.x_indices, ISMIP7_melt_rate.y_indices, grid, subdomain)
    S_loc = spatial_on_subdomain(ISMIP7_melt_rate.S_loc, grid, subdomain)
    T_loc = spatial_on_subdomain(ISMIP7_melt_rate.T_loc, grid, subdomain)
    Tf_loc = spatial_on_subdomain(ISMIP7_melt_rate.Tf_loc, grid, subdomain)

    if !(isnothing(S_loc) || isnothing(T_loc) || isnothing(Tf_loc) || isnothing(ISMIP7_melt_rate.z_forcing))
        (size(ISMIP7_melt_rate.z_forcing,1) ==
         size(S_loc,3) ==
         size(T_loc,3) ==
         size(Tf_loc,3)) ||
         throw(ArgumentError("Inconsistent size of z_forcing and ocean fields"))
    end

    return ISMIP7MeltRate(
                ISMIP7_melt_rate.so_prefix,
                ISMIP7_melt_rate.tf_prefix,
                ISMIP7_melt_rate.thetao_prefix,
                ISMIP7_melt_rate.so_varname,
                ISMIP7_melt_rate.tf_varname, 
                ISMIP7_melt_rate.thetao_varname, 
                ISMIP7_melt_rate.K,
                ISMIP7_melt_rate.shelf_slope,
                ISMIP7_melt_rate.ρ_ocean,
                ISMIP7_melt_rate.ρ_ice,
                ISMIP7_melt_rate.c_ocean,
                ISMIP7_melt_rate.L_ice,
                ISMIP7_melt_rate.β_s,
                ISMIP7_melt_rate.g,
                ISMIP7_melt_rate.f,
                S_loc,
                T_loc,
                Tf_loc,
                ISMIP7_melt_rate.z_forcing, 
                ISMIP7_melt_rate.melt_partial_cell,
                ISMIP7_melt_rate.path_to_forcing,
                xs,
                ys)
end
