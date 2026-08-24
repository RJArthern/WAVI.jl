
export NoThermoDynamics

struct NoThermoDynamics{T <: Real} <: AbstractThermoDynamics
    m_grounded :: T #uniform melt rate applied for grounded areas
end

NoThermoDynamics(; m_grounded = 0.0) = NoThermoDynamics(m_grounded)

"""
NoThermoDynamics(; <kwargs>)


Keyword arguments
=================
- m_grounded  : grounded basal melt rate (m/yr)
"""

"""
            update_ice_temperature_and_basal_melt_rate!(thermo_dynamics::NoThermoDynamics, model::AbstractModel)

no thermodynamics model is used. ice temperature is constant in time
"""

function update_ice_temperature_and_basal_melt_rate!(thermo_dynamics::NoThermoDynamics, model::AbstractModel)
    @unpack gh=model.fields
    @. gh.basal_melt = ifelse(
        gh.mask,
        gh.grounded_fraction * thermo_dynamics.m_grounded + gh.shelf_basal_melt,
        gh.basal_melt,
    )
    return model
end