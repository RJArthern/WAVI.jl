export ISMIP7SMB

using WAVI: AbstractClimateForcing
using WAVI.ClimateForcing
using NCDatasets

struct ISMIP7SMB{T <: Real,
                P <: String,
                V <: String,
                RE <: Union{AbstractArray{T,2}, Nothing},
                VSG <: Union{AbstractArray{T,2}, Nothing},
                SA <: Union{AbstractArray{T,2}, Nothing},
                RS <: Union{AbstractArray{T,2}, Nothing}} <: AbstractSurfaceMassBalance
    smb_anomaly_prefix::P
    vertical_smb_gradient_prefix::P
    smb_anomaly_varname::V 
    vertical_smb_gradient_varname::V
    reference_elevation::RE
    vertical_smb_gradient::VSG
    smb_anomaly::SA
    reference_smb::RS
    path_to_forcing::String
    # Maps local field indices to global NetCDF indices when on a subdomain.
    x_indices::Union{Nothing, UnitRange{Int}}
    y_indices::Union{Nothing, UnitRange{Int}}
end

"""
ISMIP7SMB(; <kwargs>)

Keyword arguments
=================
- climate_forcing: decriptor of climate forcing (ISMIP7_ANOMALY, ISMIP7_OCX, or ISMIP7_CONTROL)
- reference_elevation: reference elevation used for lapse rate perturbations.
- vertical_smb_gradient: gradient used for lapse rate perturbations.
- smb_anomaly: anomaly used for anomaly perturbations.    
-  reference_smb: reference smb used for anomaly perturbations. 
- path_to_forcing: path to forcing files
- x_indices: Maps local field x indices to global NetCDF x indices when on a subdomain.
- y_indices: Maps local field y indices to global NetCDF y indices when on a subdomain.
"""
function ISMIP7SMB(; 
                climate_forcing = nothing,
                reference_elevation = nothing,
                vertical_smb_gradient = nothing, 
                smb_anomaly = nothing,           
                reference_smb = nothing,
                path_to_forcing = "./",
                x_indices = nothing,
                y_indices = nothing)

    ~(climate_forcing === nothing) || throw(ArgumentError("You must pass a climate_forcing description"))

    # check that you've passed reference elevation 
    ~(reference_elevation === nothing) || throw(ArgumentError("You must pass a reference elevation file"))

    #check that you've passed a vertical smb gradient
    ~(reference_smb === nothing) || throw(ArgumentError("You must pass a reference smb"))

    return ISMIP7SMB(climate_forcing, reference_elevation,vertical_smb_gradient, smb_anomaly,reference_smb, path_to_forcing, x_indices, y_indices)

end

#Convenience Constructors
function ISMIP7SMB(climate_forcing::ISMIP7_ANOMALY,
                    reference_elevation,
                    vertical_smb_gradient, 
                    smb_anomaly,
                    reference_smb, 
                    path_to_forcing, 
                    x_indices, 
                    y_indices)

    smb_anomaly_prefix = "acabf-anomaly_"
    vertical_smb_gradient_prefix = "dacabfdz_"
    smb_anomaly_varname = "acabf-anomaly"
    vertical_smb_gradient_varname = "dacabfdz"

    return ISMIP7SMB(smb_anomaly_prefix,
                    vertical_smb_gradient_prefix,
                    smb_anomaly_varname,
                    vertical_smb_gradient_varname,
                    reference_elevation,
                    vertical_smb_gradient, 
                    smb_anomaly,
                    reference_smb, 
                    path_to_forcing, 
                    x_indices, 
                    y_indices
                    )
end

function ISMIP7SMB(climate_forcing::ISMIP7_OCX,
                    reference_elevation,
                    vertical_smb_gradient, 
                    smb_anomaly,
                    reference_smb, 
                    path_to_forcing, 
                    x_indices, 
                    y_indices)

    smb_anomaly_prefix = "acabf_"
    vertical_smb_gradient_prefix = "dacabfdz_"
    smb_anomaly_varname = "acabf"
    vertical_smb_gradient_varname = "dacabfdz"

    all(iszero.(reference_smb))||error("ISMIP7_OCX climate forcing requires zero for reference_smb")

    return ISMIP7SMB(smb_anomaly_prefix,
                    vertical_smb_gradient_prefix,
                    smb_anomaly_varname,
                    vertical_smb_gradient_varname,
                    reference_elevation,
                    vertical_smb_gradient, 
                    smb_anomaly,
                    reference_smb, 
                    path_to_forcing, 
                    x_indices, 
                    y_indices
                    )
end


function ISMIP7SMB(climate_forcing::ISMIP7_CONTROL,
                    reference_elevation,
                    vertical_smb_gradient, 
                    smb_anomaly,
                    reference_smb, 
                    path_to_forcing, 
                    x_indices, 
                    y_indices)

    smb_anomaly_prefix = "acabf_"
    vertical_smb_gradient_prefix = "dacabfdz_"
    smb_anomaly_varname = "acabf"
    vertical_smb_gradient_varname = "dacabfdz"

    all(iszero.(reference_smb))||error("ISMIP7_CONTROL climate forcing requires zero for reference_smb")

    return ISMIP7SMB(smb_anomaly_prefix,
                    vertical_smb_gradient_prefix,
                    smb_anomaly_varname,
                    vertical_smb_gradient_varname,
                    reference_elevation,
                    vertical_smb_gradient, 
                    smb_anomaly,
                    reference_smb, 
                    path_to_forcing, 
                    x_indices, 
                    y_indices
                    )
end

function reconstruct_on_grid(smb::ISMIP7SMB,grid::Grid) 
    xs, ys = grid_index_ranges(smb.x_indices, smb.y_indices, grid)
    return ISMIP7SMB(
    smb.smb_anomaly_prefix,
    smb.vertical_smb_gradient_prefix,
    smb.smb_anomaly_varname,
    smb.vertical_smb_gradient_varname,
    isnothing(smb.reference_elevation) ? zeros(grid.nx,grid.ny) : 
    size(smb.reference_elevation) == (grid.nx,grid.ny) ? smb.reference_elevation :
    throw(DimensionMismatch("Size of reference elevation is incompatible with grid")),
    isnothing(smb.vertical_smb_gradient) ? zeros(grid.nx,grid.ny) :
    size(smb.vertical_smb_gradient) == (grid.nx,grid.ny) ? smb.vertical_smb_gradient :
    throw(DimensionMismatch("Size of vertical smb gradient is incompatible with grid")),
    isnothing(smb.smb_anomaly) ? zeros(grid.nx,grid.ny) :
    size(smb.smb_anomaly) == (grid.nx,grid.ny) ? smb.smb_anomaly :
    throw(DimensionMismatch("Size of smb anomaly is incompatible with grid")),
    isnothing(smb.reference_smb) ? zeros(grid.nx,grid.ny) : 
    size(smb.reference_smb) == (grid.nx,grid.ny) ? smb.reference_smb :
    throw(DimensionMismatch("Size of reference smb is incompatible with grid")),
    smb.path_to_forcing,
    xs,
    ys)
end

function reconstruct_on_subdomain(smb::ISMIP7SMB,grid::Grid,subdomain::NTuple{4,<: Integer}) 
    xs, ys = subdomain_index_ranges(smb.x_indices, smb.y_indices, grid, subdomain)

    return ISMIP7SMB(
                    smb.smb_anomaly_prefix,
                    smb.vertical_smb_gradient_prefix,
                    smb.smb_anomaly_varname,
                    smb.vertical_smb_gradient_varname,
                    spatial_on_subdomain(smb.reference_elevation, grid, subdomain),
                    spatial_on_subdomain(smb.vertical_smb_gradient, grid, subdomain),
                    spatial_on_subdomain(smb.smb_anomaly, grid, subdomain),
                    spatial_on_subdomain(smb.reference_smb, grid, subdomain),
                    smb.path_to_forcing,
                    xs,
                    ys)
end

"""
    on_architecture(arch, smb::ISMIP7SMB)

Copy dense SMB fields onto `arch` so they can broadcast with model fields.
"""
function on_architecture(arch::AbstractArchitecture, smb::ISMIP7SMB)
    return ISMIP7SMB(
        smb.smb_anomaly_prefix,
        smb.vertical_smb_gradient_prefix,
        smb.smb_anomaly_varname,
        smb.vertical_smb_gradient_varname,
        on_architecture(arch, smb.reference_elevation),
        on_architecture(arch, smb.vertical_smb_gradient),
        on_architecture(arch, smb.smb_anomaly),
        on_architecture(arch, smb.reference_smb),
        smb.path_to_forcing,
        smb.x_indices,
        smb.y_indices,
    )
end

function update_accumulation_rate!(surface_mass_balance::ISMIP7SMB, model::AbstractModel, clock::Clock)
    @unpack reference_smb, smb_anomaly, vertical_smb_gradient, reference_elevation = surface_mass_balance
    @unpack accumulation, s = model.fields.gh

    #check the files are the right size
    (size(smb_anomaly) == (model.grid.nx, model.grid.ny)) || throw(DimensionMismatch("Size of read in smb anomaly is not compatible with grid size"))
    (size(vertical_smb_gradient) == (model.grid.nx, model.grid.ny)) || throw(DimensionMismatch("Size of read in vertical smb gradient is not compatible with grid size"))

    # Dense SMB arrays must share the architecture of `s` / `accumulation`.
    change_in_elevation = s .- reference_elevation
    accumulation .= reference_smb .+ smb_anomaly .+ vertical_smb_gradient .* change_in_elevation
    return nothing
end


function update_climate_forcing!(surface_mass_balance::ISMIP7SMB, grid::Grid, clock::Clock) 
    @unpack smb_anomaly = surface_mass_balance
    @unpack vertical_smb_gradient = surface_mass_balance
    @unpack dx = grid
    @unpack smb_anomaly_prefix, vertical_smb_gradient_prefix, path_to_forcing, x_indices, y_indices, smb_anomaly_varname, vertical_smb_gradient_varname  = surface_mass_balance

    #get the year from clock for the forcing files
    current_time = clock.time + clock.ref_time
    println(current_time)
    current_time_string = string(Int(round(current_time)))

    # load in the smb anomaly from ISMIP7
    resolution = join([string(Int(dx)), "m"])
    smb_anomaly_filename = joinpath(path_to_forcing,join([smb_anomaly_prefix,  current_time_string,".nc"]))
    smb_anomaly_ncfile   = NCDataset(smb_anomaly_filename)
    #println("read in smb anomaly forcing file: " * smb_anomaly_filename)
    @info "read in smb anomaly forcing file: $smb_anomaly_filename"


    # load in the vertical smb gradient from ISMIP7
    vertical_smb_gradient_anomaly_filename = joinpath(path_to_forcing, join([vertical_smb_gradient_prefix,  current_time_string,".nc"]))
    vertical_smb_gradient_anomaly_ncfile = NCDataset(vertical_smb_gradient_anomaly_filename)
    #println("read in vertical smb gradient forcing file: " * vertical_smb_gradient_anomaly_filename)
    @info "read in vertical smb gradient forcing file: $vertical_smb_gradient_anomaly_filename"

    is, js = forcing_index_ranges(x_indices, y_indices, smb_anomaly)
    # NetCDF reads are host arrays; copy_onto! handles host or device destinations.
    copy_onto!(
        smb_anomaly,
        replace(replace(smb_anomaly_ncfile[smb_anomaly_varname][is, js, 1], missing => NaN), NaN => 0.0),
    )
    copy_onto!(
        vertical_smb_gradient,
        replace(replace(vertical_smb_gradient_anomaly_ncfile[vertical_smb_gradient_varname][is, js, 1], missing => NaN), NaN => 0.0),
    )

    return nothing

end
