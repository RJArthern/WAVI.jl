module SlidingLaw

import WAVI.Grids: reconstruct_on_grid, reconstruct_on_subdomain
import WAVI.ClimateForcing: update_climate_forcing!
export reconstruct_on_grid, reconstruct_on_subdomain,  update_climate_forcing!
export update_β_using_sliding_law!

using Parameters

using WAVI: AbstractSlidingLaw, AbstractModel
using WAVI.Grids
using WAVI.Time


#add each of the individual sliding laws
include("./WeertmanSlidingLaw.jl")
include("./WeertmanRelaxationSlidingLaw.jl")
include("./coulomb.jl")
include("./budd.jl")
include("./tsai.jl")
include("./tsaiBudd.jl")
include("./schoof.jl")
include("./zoetIverson.jl")


"""
    update_drag_coefficient!(model::AbstractModel)

Update drag coefficient used in the sliding law to account for migration of grounding line.
This is currently only used for the Weertman sliding law and Weertman part of the Tsai sliding law, 
as they are the only sliding laws that do not directly depend on effective pressure, which already
acounts for migration of grounding line.
"""
function update_drag_coefficient!(model::AbstractModel)
    @unpack gh=model.fields
    @unpack sliding_law=model
    gh.drag_coefficient .= sliding_law.drag_coefficient .* gh.grounded_fraction
    return model
end

end