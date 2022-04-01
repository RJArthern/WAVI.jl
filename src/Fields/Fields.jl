include("HGrid.jl")
include("CGrid.jl")
include("VGrid.jl")
include("UGrid.jl")
include("SigmaGrid.jl")

"""
    Structure to hold all field variables in WAVI.jl
"""
struct Fields{T <: Real, N <: Real}
    gh::HGrid{T,N}
    gu::UGrid{T,N}
    gv::VGrid{T,N}
    gc::CGrid{T,N}
    g3d::SigmaGrid{T,N}   
    wu::UWavelets{T,N}
    wv::VWavelets{T,N}
end

"""
    setup_fields(grid, initial_conditions, solver_params, params, bed_array)

Acts as a constructor for the fields (no explicit constructor as fields `only ever called when setting up a model)
"""

function setup_fields(grid, initial_conditions, solver_params, params, bed_array)
    #Define masks for points on h-, u-, v- and c-grids that lie in model domain.
    h_mask = grid.h_mask 
    u_mask = get_u_mask(h_mask)
    v_mask = get_v_mask(h_mask)
    c_mask = get_c_mask(h_mask)

    #Remove all points on u- and v-grids with homogenous Dirichlet conditions.
    u_mask[grid.u_iszero].=false
    v_mask[grid.v_iszero].=false

    #h-grid
    #gh=HGrid(grid, params) #uncomment if using the explicit constructor method
    h =  deepcopy(initial_conditions.initial_thickness)
    ηav = deepcopy(initial_conditions.initial_viscosity[:,:,1]) #set to the viscosity on the first level for now
    gh=HGrid(
    nxh=grid.nx,
    nyh=grid.ny,
    mask=h_mask,
    b = bed_array,
    h = h,
    ηav = ηav,
    )

    #u-grid
    gu=UGrid(
    nxu=grid.nx+1,
    nyu=grid.ny,
    dx=grid.dx,
    dy=grid.dy,
    mask=u_mask,
    levels=solver_params.levels,
    u=deepcopy(initial_conditions.initial_u_veloc)
    )

    #v-grid
    gv=VGrid(
    nxv=grid.nx,
    nyv=grid.ny+1,
    dx=grid.dx,
    dy=grid.dy,
    mask=v_mask,
    levels=solver_params.levels,
    v=deepcopy(initial_conditions.initial_v_veloc)
    )

    #c-grid
    gc=CGrid(
    nxc=grid.nx-1,
    nyc=grid.ny-1,
    mask=c_mask
    )

    #3D-grid
    η = deepcopy(initial_conditions.initial_viscosity)
    θ = deepcopy(initial_conditions.initial_temperature)
    Φ = deepcopy(initial_conditions.initial_damage)
    g3d=SigmaGrid(
    nxs=grid.nx,
    nys=grid.ny,
    nσs=grid.nσ,
    σ =grid.σ,
    η = η,
    θ = θ,
    Φ = Φ,
    glen_b = glen_b.(θ,Φ,params.glen_a_ref, params.glen_n, params.glen_a_activation_energy, params.glen_temperature_ref, params.gas_const),
    quadrature_weights = grid.quadrature_weights
    )

    #Wavelet-grid, u-component.
    wu=UWavelets(nxuw=grid.nx+1,nyuw=grid.ny,levels=solver_params.levels)

    #Wavelet-grid, v-component.
    wv=VWavelets(nxvw=grid.nx,nyvw=grid.ny+1,levels=solver_params.levels)
    return Fields(gh,gu,gv,gc,g3d,wu,wv)
end