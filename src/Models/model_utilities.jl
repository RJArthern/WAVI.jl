
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