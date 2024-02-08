struct Model{T <: Real, N <: Integer,A,W, G, M <:AbstractMeltRate, PS <: AbstractParallelSpec} <: AbstractModel{T,N,M,PS}
    grid::Grid{T,N}
    params::Params{T,A,W,G}
    solver_params::SolverParams{T,N}
    initial_conditions::InitialConditions{T}
    fields::Fields{T,N}
    melt_rate::M
    parallel_spec::PS
    verbose :: Bool
end

"""
    Model(;
        grid = nothing, 
        bed_elevation = nothing,
        params = Params(),
        solver_params = SolverParams(),
        initial_conditions = InitialConditions(),
        melt_rate = UniformMeltRate(),
        parallel_spec = BasicParallelSpec(),
        verbose = true)

Construct a WAVI.jl model object.

Keyword arguments
=================

    - `grid`: (required) an instance of a `Grid` object, which defines the computational grid
    - `bed_elevation`: (required) an array of size `grid.nx` x `grid.ny` which defines the bed elevation
    - `params`: a `Params` object that defines physical parameters 
    - `solver_params`: a `SolverParams` object that defines parameters relating to the numerical scheme
    - `initial_conditions`: an `InitialConditions` object that (optionally) defines the initial ice thickness, temperature, viscosity, and damage
    - `melt_rate`: a melt rate model, responsible for controlling/setting the basal melt rate
    - `parallel_spec`: specification of parallel computation method.
    - `verbose`: specifies whether to output information about the solve and residuals

"""
function Model(;
    grid = nothing, 
    bed_elevation = nothing,
    params = Params(),
    solver_params = SolverParams(),
    initial_conditions = InitialConditions(),
    melt_rate = UniformMeltRate(),
    parallel_spec = BasicParallelSpec(),
    verbose = true)

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

    ## Parameter fields checks 
    #if weertman c passed as a scalar, replace weertman_c parameters with matrix of this value
    if isa(params.weertman_c, Number) 
        params = @set params.weertman_c = params.weertman_c*ones(grid.nx,grid.ny)
    end
    #check size compatibility of resulting weertman C
    (size(params.weertman_c)==(grid.nx,grid.ny)) || throw(DimensionMismatch("Size of input weertman c must match grid size (i.e. $(grid.nx) x $(grid.ny))"))
    
    #if accumulation is passed as a scalar, replace accumulation parameters with matrix of this value
    if isa(params.accumulation_rate, Number) 
        params = @set params.accumulation_rate = params.accumulation_rate*ones(grid.nx,grid.ny)
    end
    #check size compatibility of resulting accumulation rate
    (size(params.accumulation_rate)==(grid.nx,grid.ny)) || throw(DimensionMismatch("Size of input accumulation must match grid size (i.e. $(grid.nx) x $(grid.ny))"))

    #if accumulation is passed as a scalar, replace accumulation parameters with matrix of this value
    if isa(params.glen_a_ref, Number) 
        params = @set params.glen_a_ref = params.glen_a_ref*ones(grid.nx,grid.ny)
    end
    #check size compatibility of resulting glen a ref
    (size(params.glen_a_ref)==(grid.nx,grid.ny)) || throw(DimensionMismatch("Size of input glen_a_ref must match grid size (i.e. $(grid.nx) x $(grid.ny))"))




    #Setup the fields 
    fields = setup_fields(grid, initial_conditions, solver_params, params, bed_array)

    #Use type constructor to build initial state with no extra physics
    model=Model(grid,params,solver_params,initial_conditions,fields,melt_rate,parallel_spec, verbose)

    return model
end

include("model_utilities.jl")
include("update_state.jl")
include("update_velocities.jl")