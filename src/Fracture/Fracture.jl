module Fracture

import WAVI.Grids: reconstruct_on_grid, reconstruct_on_subdomain
import WAVI.ClimateForcing: update_climate_forcing!
export reconstruct_on_grid, reconstruct_on_subdomain,  update_climate_forcing!
export update_damage!, update_strain_history!

using Parameters, NCDatasets

using WAVI: AbstractFracture, AbstractModel
using WAVI.Advection
using WAVI.Time
using WAVI.Grids
using WAVI.Architectures: AbstractArchitecture
import WAVI.Architectures: on_architecture

# Default: most fracture schemes have no dense arrays to move.
on_architecture(::AbstractArchitecture, fracture::AbstractFracture) = fracture

get_fracture(model::AbstractModel{T,N}) where {T,N} = model.fracture

"""
    update_damage!(model::AbstractModel)

Update the damage field and the strain history from 'ice shelf' parts of strain rate tensor, neglecting all vertical shear.
"""
update_damage!(model::AbstractModel{T,N};kwargs ...) where {T,N} = update_damage!(get_fracture(model),model; kwargs ...)


"""
    update_strain_history!(model::AbstractModel)

Update the strain history stored in the model structure.
"""
update_strain_history!(model::AbstractModel{T,N};kwargs ...) where {T,N} = update_strain_history!(get_fracture(model),model; kwargs ...)

include("./ConstantDamage.jl")
include("./DruckerPragerPhaseField.jl")
include("./ISMIP7hydrofracture.jl")
include("./utils.jl")

end