module SurfaceMassBalance

import WAVI.Grids: reconstruct_on_grid, reconstruct_on_subdomain
import WAVI.ClimateForcing: update_climate_forcing!
export reconstruct_on_grid, reconstruct_on_subdomain,  update_climate_forcing!
export update_accumulation_rate!

using WAVI: AbstractSurfaceMassBalance, AbstractModel, AbstractClimateForcing
using Parameters
using WAVI.Time
using WAVI.Grids
using WAVI.Utilities: copy_onto!
using WAVI.Architectures: AbstractArchitecture
import WAVI.Architectures: on_architecture

# Default: most SMB schemes have no dense arrays to move.
on_architecture(::AbstractArchitecture, smb::AbstractSurfaceMassBalance) = smb

#return the sceme used to compute surface mass balance
get_surface_mass_balance(model::AbstractModel) = model.surface_mass_balance

#dispatch to specialised methods depending on surface mass balance law
"""

    update_accumulation_rate!(model::AbstractModel, clock::Clock)

Update the accumulation rate in the model using the time from clock.
"""
update_accumulation_rate!(model::AbstractModel, clock::Clock) = update_accumulation_rate!(get_surface_mass_balance(model)::AbstractSurfaceMassBalance, model, clock)

#include all specialised methods for computing accumulation rate
include("./AccumulationFromParams.jl")
include("./UniformAccumulationRate.jl")
include("./BinfileAccumulationRate.jl")
include("./ISMIP7SMB.jl")

#default behaviour if no clock passed
update_accumulation_rate!(model::AbstractModel) = update_accumulation_rate!(model::AbstractModel, Clock()) 

end