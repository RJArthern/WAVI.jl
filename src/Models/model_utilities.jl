
"""
check_initial_conditions(initial_conditions, params)

Check whether initial conditions have been specified. Default them to standard values if not
"""
function check_initial_conditions(initial_conditions, params, grid)
if all(isnan.(initial_conditions.initial_thickness))
    default_thickness = params.default_thickness
    #@info "Did not find a specified initial thickness, reverting to default value specified in params ($default_thickness m everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
    initial_conditions = @set initial_conditions.initial_thickness =  default_thickness*ones(grid.nx, grid.ny)
end

if all(isnan.(initial_conditions.initial_grounded_fraction))
    #@info "Did not find a specified grounded fraction, reverting to default value of one everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
    initial_conditions = @set initial_conditions.initial_grounded_fraction =  ones(grid.nx, grid.ny)
end

if all(isnan.(initial_conditions.initial_u_veloc))
    #@info "Did not find a specified initial u velocity, reverting to default value of zero everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
    initial_conditions = @set initial_conditions.initial_u_veloc =  zeros(grid.nx+1, grid.ny)
end

if all(isnan.(initial_conditions.initial_v_veloc))
    #@info "Did not find a specified initial u velocity, reverting to default value of zero everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
    initial_conditions = @set initial_conditions.initial_v_veloc =  zeros(grid.nx, grid.ny+1)
end

if all(isnan.(initial_conditions.initial_viscosity))
    default_viscosity = params.default_viscosity
    #@info "Did not find a specified initial viscosity, reverting to default value specified in params ($default_viscosity Pa s everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
    initial_conditions = @set initial_conditions.initial_viscosity =  default_viscosity*ones(grid.nx, grid.ny, grid.nσ)
end

if all(isnan.(initial_conditions.initial_temperature))
    default_temperature = params.default_temperature
    #@info "Did not find a specified initial temperature, reverting to default value specified in params ($default_temperature K everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
    initial_conditions = @set initial_conditions.initial_temperature =  default_temperature*ones(grid.nx, grid.ny, grid.nσ)
end

if all(isnan.(initial_conditions.initial_damage))
    default_damage = params.default_damage
    #@info "Did not find a specified initial damage field, reverting to default value specified in params ($default_damage everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
    initial_conditions = @set initial_conditions.initial_damage =  default_damage*ones(grid.nx, grid.ny, grid.nσ)
end



#check sizes are compatible
(size(initial_conditions.initial_thickness) == (grid.nx, grid.ny)) || throw(DimensionMismatch("Initial thickness field is not compatible with grid size. Input thickess field is has size $(size(initial_conditions.initial_thickness)), which must match horizontal grid size ($(grid.nx) x $(grid.ny))"))
(size(initial_conditions.initial_grounded_fraction) == (grid.nx, grid.ny)) || throw(DimensionMismatch("Initial grounded fraction field is not compatible with grid size. Input grounded fraction field is has size $(size(initial_conditions.initial_grounded_fraction)), which must match horizontal grid size ($(grid.nx) x $(grid.ny))"))
(size(initial_conditions.initial_u_veloc) == (grid.nx+1, grid.ny)) || throw(DimensionMismatch("Initial u_velocity field is not compatible with grid size. Input u_velocity field has size $(size(initial_conditions.initial_u_veloc)), which must match horizontal grid size ($(grid.nx+1) x $(grid.ny))"))
(size(initial_conditions.initial_v_veloc) == (grid.nx, grid.ny+1)) || throw(DimensionMismatch("Initial v_velocity field is not compatible with grid size. Input v_velocity field has size $(size(initial_conditions.initial_v_veloc)), which must match horizontal grid size ($(grid.nx) x $(grid.ny+1))"))
(size(initial_conditions.initial_temperature) == (grid.nx, grid.ny, grid.nσ)) || throw(DimensionMismatch("Initial temperature field is not compatible with grid size. Input temperature field is has size $(size(initial_conditions.initial_temperature)), which must match 3D grid size ($(grid.nx), $(grid.ny), $(grid.nσ))"))
(size(initial_conditions.initial_viscosity) == (grid.nx, grid.ny, grid.nσ)) || throw(DimensionMismatch("Initial viscosity field is not compatible with grid size. Input temperature field is has size $(size(initial_conditions.initial_temperature)), which must match 3D grid size ($(grid.nx), $(grid.ny), $(grid.nσ))"))
(size(initial_conditions.initial_damage) == (grid.nx, grid.ny, grid.nσ)) || throw(DimensionMismatch("Initial damage field is not compatible with grid size.Input temperature field is has size $(size(initial_conditions.initial_temperature)), which must match 3D grid size ($(grid.nx), $(grid.ny), $(grid.nσ))"))


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
