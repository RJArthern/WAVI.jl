#Struct to hold model state
include("HGrid.jl")
include("CGrid.jl")
include("VGrid.jl")
include("UGrid.jl")
include("SigmaGrid.jl")

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
    ηav = deepcopy(initial_conditions.initial_viscosity)
    gh=HGrid(
    Nx=grid.nx,
    Ny=grid.ny,
    mask=h_mask,
    b = bed_array,
    h = h,
    ηav = ηav,
    )

    #u-grid
    gu=UGrid(
    Nx=grid.nx+1,
    Ny=grid.ny,
    dx=grid.dx,
    dy=grid.dy,
    mask=u_mask,
    levels=solver_params.levels
    )

    #v-grid
    gv=VGrid(
    Nx=grid.nx,
    Ny=grid.ny+1,
    dx=grid.dx,
    dy=grid.dy,
    mask=v_mask,
    levels=solver_params.levels
    )

    #c-grid
    gc=CGrid(
    x0=grid.x0,
    y0=grid.y0,
    nx=grid.nx-1,
    ny=grid.ny-1,
    dx=grid.dx,
    dy=grid.dy,
    mask=c_mask
    )

    #3D-grid
    g3d=SigmaGrid(
    nx=grid.nx,
    ny=grid.ny,
    nσ=grid.nσ,
    η = fill(params.default_viscosity,grid.nx,grid.ny,grid.nσ),
    θ = fill(params.default_temperature,grid.nx,grid.ny,grid.nσ),
    Φ = fill(params.default_damage,grid.nx,grid.ny,grid.nσ),
    glen_b = fill(glen_b(params.default_temperature,params.default_damage,params),grid.nx,grid.ny,grid.nσ)
    )

    #Wavelet-grid, u-component.
    wu=UWavelets(nx=grid.nx+1,ny=grid.ny,levels=solver_params.levels)

    #Wavelet-grid, v-component.
    wv=VWavelets(nx=grid.nx,ny=grid.ny+1,levels=solver_params.levels)
    return Fields(gh,gu,gv,gc,g3d,wu,wv)
end
