struct Model{T <: Real, N <: Integer} <: AbstractModel{T,N}
    grid::Grid{T,N}
    params::Params{T}
    solver_params::SolverParams{T,N}
    initial_conditions::InitialConditions{T}
    fields::Fields{T,N}
    extra_physics::Dict{String, Any}
end

"""
 Model constructor
"""
function Model(;
    grid = nothing, 
    bed_elevation = nothing,
    params = Params(),
    solver_params = SolverParams(),
    initial_conditions = InitialConditions(),
    extra_physics = Dict{String, Any}())

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
    #check size compatibility of resulting weertman C
    (size(params.accumulation_rate)==(grid.nx,grid.ny)) || throw(DimensionMismatch("Size of input weertman c must match grid size (i.e. $(grid.nx) x $(grid.ny))"))

    #Setup the fields 
    fields = setup_fields(grid, initial_conditions, solver_params, params, bed_array)
    
    #Use type constructor to build initial state with no extra physics
    model=Model(grid,params,solver_params,initial_conditions,fields,extra_physics)

    return model
end

include("model_utilities.jl")
include("update_state.jl")
include("update_velocities.jl")