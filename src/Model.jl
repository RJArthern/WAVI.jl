
#Struct to hold model state
struct Model{T <: Real, N <: Integer} <: AbstractModel{T,N}
    grid::Grid{T,N}
    params::Params{T}
    solver_params::SolverParams{T,N}
    initial_conditions::InitialConditions{T}
    gh::HGrid{T,N}
    gu::UGrid{T,N}
    gv::VGrid{T,N}
    gc::CGrid{T,N}
    g3d::SigmaGrid{T,N}   
    wu::UWavelets{T,N}
    wv::VWavelets{T,N}
end

"""
 Model constructor
"""
function Model(;
    grid = nothing, 
    bed_elevation = nothing,
    params = Params(),
    solver_params = SolverParams(),
    initial_conditions = InitialConditions())

    #check that a grid and bed has been inputted
    ~(grid === nothing) || throw(ArgumentError("You must specify an input grid"))
    ~(bed_elevation === nothing) || throw(ArgumentError("You must input a bed elevation"))
    
    #check types
    #if a functional bed has been specified, convert to an array
    bed_array = zeros(grid.nx,grid.nx) #initialize a bed array
    try
    bed_array = get_bed_elevation(bed_elevation, grid)
    catch
    throw(ArgumentError("bed elevation must be of type function or array"))
    end
            #check the size of the bed
    #(Base.size(bed_array) = (grid.nx, grid.ny)) || throw(ArgumentError("Bed and grid sizes must be identical"))
    

    #Check initial conditions, and revert to default values if not
    initial_conditions = check_initial_conditions(initial_conditions, params, grid)

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
    #h =  deepcopy(initial_conditions.initial_thickness)
    h =  deepcopy(initial_conditions.initial_thickness)
    ηav = deepcopy(initial_conditions.initial_viscosity)
    gh=HGrid(
    x0=grid.x0,
    y0=grid.y0,
    nx=grid.nx,
    ny=grid.ny,
    dx=grid.dx,
    dy=grid.dy,
    mask=h_mask,
    b = bed_array,
    h = h,
    ηav = ηav,
    )

    #u-grid
    gu=UGrid(
    x0=grid.x0,
    y0=grid.y0,
    nx=grid.nx+1,
    ny=grid.ny,
    dx=grid.dx,
    dy=grid.dy,
    mask=u_mask,
    levels=solver_params.levels
    )

    #v-grid
    gv=VGrid(
    x0=grid.x0,
    y0=grid.y0,
    nx=grid.nx,
    ny=grid.ny+1,
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

    #Use type constructor to build initial state.
    wavi=Model(grid,params,solver_params,initial_conditions,gh,gu,gv,gc,g3d,wu,wv)

    return wavi
end


"""
    check_initial_conditions(initial_conditions, params)

Check whether initial conditions have been specified. Default them to standard values if not
"""
function check_initial_conditions(initial_conditions, params, grid)
    if initial_conditions.initial_thickness == zeros(10,10)
        default_thickness = params.default_thickness
        println("Did not find a specified initial thickness, reverting to default value specified in params ($default_thickness m everywhere)")
        initial_conditions = @set initial_conditions.initial_thickness =  default_thickness*ones(grid.nx, grid.ny)
    end
    
    if initial_conditions.initial_viscosity == zeros(10,10)
        default_viscosity = params.default_viscosity
        println("Did not find a specified initial viscosity, reverting to default value specified in params ($default_viscosity Pa s everywhere)")
        initial_conditions = @set initial_conditions.initial_viscosity =  default_viscosity*ones(grid.nx, grid.ny)
    end

    if initial_conditions.initial_temperature == zeros(10,10)
        default_temperature = params.default_temperature
        println("Did not find a specified initial temperature, reverting to default value specified in params ($default_temperature K everywhere)")
        initial_conditions = @set initial_conditions.initial_temperature =  default_temperature*ones(grid.nx, grid.ny)
    end

    if initial_conditions.initial_damage == zeros(10,10)
        default_damage = params.default_damage
        println("Did not find a specified initial damage field, reverting to default value specified in params ($default_damage everywhere)")
        initial_conditions = @set initial_conditions.initial_damage =  default_damage*ones(grid.nx, grid.ny)
    end
    return initial_conditions
end

# Utility functions
function get_bed_elevation(bed_elevation::F, grid) where (F <: Function)
    bed_array = bed_elevation.(grid.xxh, grid.yyh)
    return bed_array
end

function get_bed_elevation(bed_elevation::Array{T,2}, grid) where (T <: Real)
    bed_array = bed_elevation
    return bed_array
end
