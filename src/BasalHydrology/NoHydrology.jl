export NoHydrology

struct NoHydrology <: AbstractBasalHydrology end

"""
            update_basal_water_thickness_effective_pressure!(basal_hydrology::NoHydrology, model::AbstractModel; update_basal_water_thickness::Bool = true)

no basal hydrology model is used. effective pressure is constant in time
"""

function update_basal_water_thickness_effective_pressure!(basal_hydrology::NoHydrology, model::AbstractModel; update_basal_water_thickness::Bool = true)
    @unpack gh=model.fields
    copy_onto!(gh.effective_pressure, model.params.effective_pressure)
    gh.effective_pressure .*= gh.grounded_fraction
    return model
end