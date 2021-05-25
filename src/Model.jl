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

    #Setup the fields 
    fields = setup_fields(grid, initial_conditions, solver_params, params, bed_array)
    
    #Use type constructor to build initial state.
    model=Model(grid,params,solver_params,initial_conditions,fields,extra_physics)

    return model
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


