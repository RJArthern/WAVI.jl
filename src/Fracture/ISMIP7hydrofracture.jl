export ISMIP7Hydrofracture

using WAVI: AbstractClimateForcing
using WAVI.ClimateForcing
using WAVI.Utilities: copy_onto!
using WAVI.Architectures: AbstractArchitecture
import WAVI.Architectures: on_architecture
using NCDatasets

struct ISMIP7Hydrofracture{P <: String, V<:String, T <: Real, CM <: Union{AbstractArray{T,2}, Nothing}} <: AbstractFracture
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

  println("read in fracture mask file: " * ice_shelf_collapse_mask_filename)
  @info "read in fracture mask file: $ice_shelf_collapse_mask_filename"

  if isnothing(ice_shelf_collapse_mask)
      throw(ArgumentError("ISMIP7 hydrofracture mask has not been allocated; reconstruct on the grid first"))
  end
  is, js = forcing_index_ranges(x_indices, y_indices, ice_shelf_collapse_mask)
  # NetCDF mask is typically Int8; destination is Float64 (host or device).
  mask_host = Float64.(ice_shelf_collapse_mask_ncfile[hydrofracture_varname][is, js, 1])
  copy_onto!(ice_shelf_collapse_mask, mask_host)

  return nothing
end

"""
    on_architecture(arch, fracture::ISMIP7Hydrofracture)

Copy the collapse mask onto `arch` so it can broadcast with model fields.
"""
function on_architecture(arch::AbstractArchitecture, fracture::ISMIP7Hydrofracture)
    return ISMIP7Hydrofracture(
        fracture.hydrofracture_prefix,
        fracture.hydrofracture_varname,
        fracture.damage_value,
        fracture.partially_floating_cells,
        on_architecture(arch, fracture.ice_shelf_collapse_mask),
        fracture.path_to_forcing,
        fracture.x_indices,
        fracture.y_indices,
    )
end

function update_damage!(fracture::ISMIP7Hydrofracture,model::AbstractModel{T,N};kwargs...) where {T,N}
  @unpack gh,g3d = model.fields

  if fracture.partially_floating_cells
    g3d.Φ .= max.(g3d.Φ,fracture.ice_shelf_collapse_mask .* (1 .- gh.grounded_fraction) .* fracture.damage_value)
  else
    g3d.Φ .= max.(g3d.Φ,fracture.ice_shelf_collapse_mask .* (gh.grounded_fraction .== 0.0) .* fracture.damage_value)
  end

  return model
end

function update_strain_history!(fracture::ISMIP7Hydrofracture,model::AbstractModel{T,N};kwargs...) where {T,N}
  return model
end

function reconstruct_on_grid(fracture::ISMIP7Hydrofracture, grid::Grid)
    xs, ys = grid_index_ranges(fracture.x_indices, fracture.y_indices, grid)
    return ISMIP7Hydrofracture(
        fracture.hydrofracture_prefix,
        fracture.hydrofracture_varname,
        fracture.damage_value,
        fracture.partially_floating_cells,
        isnothing(fracture.ice_shelf_collapse_mask) ? zeros(grid.nx,grid.ny) :
        size(fracture.ice_shelf_collapse_mask) == (grid.nx,grid.ny) ? fracture.ice_shelf_collapse_mask :
        throw(DimensionMismatch("Size of ice shelf collapse_mask is incompatible with grid")),
        fracture.path_to_forcing,
        xs,
        ys)
end

function reconstruct_on_subdomain(fracture::ISMIP7Hydrofracture, grid::Grid,subdomain::NTuple{4,<: Integer})
    xs, ys = subdomain_index_ranges(fracture.x_indices, fracture.y_indices, grid, subdomain)
    return ISMIP7Hydrofracture(
        fracture.hydrofracture_prefix,
        fracture.hydrofracture_varname,
        fracture.damage_value,
        fracture.partially_floating_cells,
        spatial_on_subdomain(fracture.ice_shelf_collapse_mask, grid, subdomain),
        fracture.path_to_forcing,
        xs,
        ys)
end
