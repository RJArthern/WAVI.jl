export ISMIP7Hydrofracture

using WAVI: AbstractClimateForcing
using WAVI.ClimateForcing
using NCDatasets

struct ISMIP7Hydrofracture{P <: String, V<:String, T <: Real, CM <: Union{Array{T,2}, Nothing}} <: AbstractFracture
      hydrofracture_prefix::P
      hydrofracture_varname::V
      damage_value :: T
      partially_floating_cells :: Bool
      ice_shelf_collapse_mask :: CM
      path_to_forcing :: String
      # Maps local mask indices to global NetCDF indices when on a subdomain.
      x_indices::Union{Nothing, UnitRange{Int}}
      y_indices::Union{Nothing, UnitRange{Int}}
end

"""
ISMIP7Hydrofracture(; <kwargs>)


Keyword arguments
=================
- climate_forcing           : decriptor of climate forcing (ISMIP7_ANOMALY, ISMIP7_OCX, or ISMIP7_CONTROL)
- damage_value              : damage value of all floating grid cells within the ice shelf collapse mask
- partially_floating_cells  : consider partially floating cells?
- ice_shelf_collapse_mask   : mask is 0 when ice shelves are sustainable and 1 when they collapse because excessive amounts of meltwater
- path_to_forcing           : path to the forcing files folder
- x_indices                 : Maps local field x indices to global NetCDF x indices when on a subdomain.
- y_indices                 : Maps local field y indices to global NetCDF y indices when on a subdomain.
"""

function ISMIP7Hydrofracture(;
              climate_forcing = nothing,
              damage_value = 0.99,
              partially_floating_cells = false,
              ice_shelf_collapse_mask = nothing,
              path_to_forcing = "./",
              x_indices = nothing,
              y_indices = nothing
              )

        
    ~(climate_forcing === nothing) || throw(ArgumentError("You must pass a climate_forcing description"))

    return ISMIP7Hydrofracture(climate_forcing,
              damage_value,
              partially_floating_cells,
              ice_shelf_collapse_mask, 
              path_to_forcing,
              x_indices,
              y_indices)
end

#Convenience Constructors
function ISMIP7Hydrofracture(climate_forcing::T,
              damage_value,
              partially_floating_cells,
              ice_shelf_collapse_mask, 
              path_to_forcing,
              x_indices,
              y_indices) where {T <: Union{ISMIP7_ANOMALY,ISMIP7_OCX,ISMIP7_CONTROL}}

    hydrofracture_prefix = "ice_shelf_collapse_mask_"
    hydrofracture_varname = "mask"

    return ISMIP7Hydrofracture(hydrofracture_prefix,
                    hydrofracture_varname,
                    damage_value,
                    partially_floating_cells,
                    ice_shelf_collapse_mask, 
                    path_to_forcing,
                    x_indices,
                    y_indices)
end

function update_climate_forcing!(fracture::ISMIP7Hydrofracture, grid::Grid, clock::Clock)
  @unpack dx = grid
  @unpack hydrofracture_prefix, path_to_forcing, ice_shelf_collapse_mask, x_indices, y_indices,hydrofracture_varname = fracture


  #get the year from clock for the forcing files
  current_time = clock.time + clock.ref_time
  current_time_string = string(Int(round(current_time)))

  # load in the ice shelf collapse mask from ISMIP7
  resolution = join([string(Int(dx)), "m"])
  ice_shelf_collapse_mask_filename = joinpath(path_to_forcing, join([hydrofracture_prefix,  current_time_string,".nc"]))
  ice_shelf_collapse_mask_ncfile   = NCDataset(ice_shelf_collapse_mask_filename)
  mask_full = ice_shelf_collapse_mask_ncfile[hydrofracture_varname][:,:,1]

  println("read in fracture mask file: " * ice_shelf_collapse_mask_filename)
  @info "read in fracture mask file: $ice_shelf_collapse_mask_filename"

  if size(ice_shelf_collapse_mask) == size(mask_full)
      ice_shelf_collapse_mask .= mask_full
  else
      is = isnothing(x_indices) ? (1:size(ice_shelf_collapse_mask, 1)) : x_indices
      js = isnothing(y_indices) ? (1:size(ice_shelf_collapse_mask, 2)) : y_indices
      ice_shelf_collapse_mask .= mask_full[is, js]
  end

  return nothing
end

function update_damage!(fracture::ISMIP7Hydrofracture,model::AbstractModel{T,N};kwargs...) where {T,N}
  @unpack gh,g3d = model.fields

  if fracture.partially_floating_cells
    change_mask = Bool.(fracture.ice_shelf_collapse_mask)
    g3d.Φ[findall(change_mask),:] .= g3d.Φ[findall(change_mask),:].*gh.grounded_fraction[findall(change_mask)] .+
        (one(eltype(gh.grounded_fraction)) .- gh.grounded_fraction[findall(change_mask)]) .* 
                  max.(g3d.Φ[findall(change_mask),:], fracture.damage_value)
  else
    change_mask = Bool.(fracture.ice_shelf_collapse_mask) .& (gh.grounded_fraction .== 0.0)
    g3d.Φ[findall(change_mask),:] .= max.(g3d.Φ[findall(change_mask),:],fracture.damage_value)
  end

  return model
end

function update_strain_history!(fracture::ISMIP7Hydrofracture,model::AbstractModel{T,N};kwargs...) where {T,N}
  return model
end

function reconstruct_on_grid(fracture::ISMIP7Hydrofracture, grid::Grid)
    return ISMIP7Hydrofracture(
        fracture.hydrofracture_prefix,
        fracture.hydrofracture_varname,
        fracture.damage_value,
        fracture.partially_floating_cells,
        isnothing(fracture.ice_shelf_collapse_mask) ? zeros(grid.nx,grid.ny) :
        size(fracture.ice_shelf_collapse_mask) == (grid.nx,grid.ny) ? fracture.ice_shelf_collapse_mask :
        throw(DimensionMismatch("Size of ice shelf collapse_mask is incompatible with grid")),
        fracture.path_to_forcing,
        1:grid.nx,
        1:grid.ny)
end

function reconstruct_on_subdomain(fracture::ISMIP7Hydrofracture, grid::Grid,subdomain::NTuple{4,<: Integer})
    x_start, x_end, y_start, y_end = subdomain
    parent_x = isnothing(fracture.x_indices) ? (1:size(fracture.ice_shelf_collapse_mask, 1)) : fracture.x_indices
    parent_y = isnothing(fracture.y_indices) ? (1:size(fracture.ice_shelf_collapse_mask, 2)) : fracture.y_indices
    xs = parent_x[x_start:x_end]
    ys = parent_y[y_start:y_end]
    return ISMIP7Hydrofracture(
        fracture.hydrofracture_prefix,
        fracture.hydrofracture_varname,
        fracture.damage_value,
        fracture.partially_floating_cells,
        size(fracture.ice_shelf_collapse_mask) == size(grid)[1:2] ? fracture.ice_shelf_collapse_mask[x_start:x_end, y_start:y_end] : fracture.ice_shelf_collapse_mask,
        fracture.path_to_forcing,
        xs,
        ys)
end
