struct UniformMeltFloatOnlyBasinSpecific{T <: Real} <: AbstractMeltRate 
    melt_constant_basin_21  :: T
    melt_constant_basin_22  :: T
    melt_partial_cell :: Bool
    ρi :: T
    ρw :: T
end

"""
    function UniformMeltFloatOnlyBasinSpecific(; <kwargs>)

Construct a melt rate to specify the melt rate as a constant under floating ice and zero elsewhere

Keyword arguments
=================
- melt_constant_basin_21 : uniform melt rate in basin 21 (Thwaites) (defaults to zero)
- melt_constant_basin_22 : uniform melt rate in basin 22 (PIG) (defaults to zero)
- melt_partial_cell : whether to melt under partially grounded cells or not (default false)
- ρi : Ice density
- ρw : Water density
"""
UniformMeltFloatOnlyBasinSpecific(; melt_constant_basin_21 = 0.0, melt_constant_basin_22 = 0.0, melt_partial_cell= false, ρi = 918.0, ρw = 1028.0) = UniformMeltFloatOnlyBasinSpecific(melt_constant_basin_21,melt_constant_basin_22, melt_partial_cell,ρi, ρw)

"""
    update_melt_rate(melt_rate::UniformMeltFloatOnly, fields, grid) 

Update the melt rate when for the UniformMeltFloatOnlyBasinSpecific type
"""
function update_melt_rate!(melt_rate::UniformMeltFloatOnly, fields, grid) 
  #  @unpack basal_melt, h, b, grounded_fraction = fields.gh #get the ice thickness and grounded fraction

  m = zeros(grid.nx,grid.ny);

    if (melt_rate.melt_partial_cell)  #partial cell melting 
        m[(basin_ID .==21)] .=   melt_rate.melt_constant_basin_21.* (1 .- fields.gh.grounded_fraction[(basin_ID .==21)])
        m[(basin_ID .==22)] .=   melt_rate.melt_constant_basin_22.* (1 .- fields.gh.grounded_fraction[(basin_ID .==22)])
    elseif ~melt_rate.melt_partial_cell #no partial cell melting
        m[(basin_ID .==21)] .=  melt_rate.melt_constant_basin_21.* (1 .- fields.gh.grounded_fraction[(basin_ID .==21)])
        m[(basin_ID .==22)] .=  melt_rate.melt_constant_basin_22.* (1 .- fields.gh.grounded_fraction[(basin_ID .==22)])
        m[.~(fields.gh.grounded_fraction .== 0)] .= 0 
    end
    
    fields.gh.basal_melt[:] .= m[:]
end

