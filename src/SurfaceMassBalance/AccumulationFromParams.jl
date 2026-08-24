export AccumulationFromParams

struct AccumulationFromParams <: AbstractSurfaceMassBalance end


function update_accumulation_rate!(surface_mass_balance::AccumulationFromParams, model::AbstractModel, clock::Clock)
    @unpack params = model
    @unpack gh=model.fields
    copy_onto!(gh.accumulation, params.accumulation_rate)
    return model
end

