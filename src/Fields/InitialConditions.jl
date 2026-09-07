using Parameters

"""
InitialConditions(; 
                    initial_thickness = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_grounded_fraction = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_u_veloc = fill!(Array{Float64}(undef,1,1),NaN)
                    initial_v_veloc = fill!(Array{Float64}(undef,1,1),NaN)
                    initial_viscosity = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_temperature = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_damage = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_strain_history = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_u_veloc=fill!(Array{Float64}(undef,1,1),NaN),
                    initial_v_veloc=fill!(Array{Float64}(undef,1,1),NaN),
                    initial_basal_water_thickness = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_hydraulic_potential_b = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_effective_pressure = fill!(Array{Float64}(undef,1,1),NaN),
                    initial_θ_ave = fill!(Array{Float64}(undef,1,1),NaN)
                    initial_preBfactor = fill!(Array{Float64}(undef,1,1),NaN))

Construct a WAVI.jl initial conditions object. 
Unpassed arguments default to 1x1 nan matrix; unspecified initial conditions are overwritten by default values specified in Params structure when model is assembled.

Keyword arguments
=================

- 'initial_thickness': (nx x ny) matrix defining ice thickness at t = 0
- 'initial_grounded_fraction': (nx x ny) matrix defining grounded fraction at t = 0
- 'initial_u_veloc' : Initial guess for depth-averaged u velocity, including any fixed velocities (nx+1,ny). 
- 'initial_v_veloc' : Initial guess for depth-averaged v velocity, including any fixed velocities (nx,ny+1).
- 'initial_viscosity': (nx x ny x nz) matrix defining viscosity on sigma levels at t = 0
- 'initial_temperature': (nx x ny x nz) matrix defining temperature on sigma levels at t = 0
- 'initial_damage': (nx x ny x nz) matrix defining ice damage at t = 0
- 'initial_strain_history': (nx x ny x nz) matrix defining maximum strain energy previously encountered at t = 0
- 'initial_basal_water_thickness': (nx x ny) matrix defining basal water thickness at t = 0
- 'initial_hydraulic_potential_b' : (nx x ny) matrix defining hydraulic potential at bed at t = 0
- 'initial_effective_pressure' : (nx x ny) matrix defining effective pressure at t = 0
- 'initial_θ_ave' : (nx x ny) matrix defining depth-averaged temperature at t = 0
- 'initial_preBfactor' : (nx x ny) matrix defining preBfactor at t = 0
"""
# TODO: this should be in Configuration.jl? 
@with_kw struct InitialConditions{T <: Real}
    initial_thickness::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_grounded_fraction::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_u_veloc::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_v_veloc::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_viscosity::Array{T,3} = fill!(Array{Float64}(undef,1,1,1),NaN)
    initial_temperature::Array{T,3} = fill!(Array{Float64}(undef,1,1,1),NaN)
    initial_damage::Array{T,3} = fill!(Array{Float64}(undef,1,1,1),NaN)
    initial_strain_history::Array{T,3} = fill!(Array{Float64}(undef,1,1,1),NaN)
    initial_basal_water_thickness::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_hydraulic_potential_b::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_effective_pressure::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_θ_ave::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
    initial_preBfactor::Array{T,2} = fill!(Array{Float64}(undef,1,1),NaN)
end


function reconstruct_on_grid(initial_conditions::InitialConditions, params::Params, grid::Grid)


    if all(isnan.(initial_conditions.initial_thickness))
        default_thickness = params.default_thickness
        #@info "Did not find a specified initial thickness, reverting to default value specified in params ($default_thickness m everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_thickness =  default_thickness*ones(grid.nx, grid.ny)
    else
        initial_thickness = initial_conditions.initial_thickness
    end

    if all(isnan.(initial_conditions.initial_grounded_fraction))
        #@info "Did not find a specified grounded fraction, reverting to default value of one everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_grounded_fraction =  ones(grid.nx, grid.ny)
    else
        initial_grounded_fraction = initial_conditions.initial_grounded_fraction
    end

    if all(isnan.(initial_conditions.initial_u_veloc))
        #@info "Did not find a specified initial u velocity, reverting to default value of zero everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_u_veloc =  zeros(grid.nx+1, grid.ny)
    else
        initial_u_veloc = initial_conditions.initial_u_veloc
    end

    if all(isnan.(initial_conditions.initial_v_veloc))
        #@info "Did not find a specified initial u velocity, reverting to default value of zero everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_v_veloc =  zeros(grid.nx, grid.ny+1)
    else
        initial_v_veloc = initial_conditions.initial_v_veloc
    end

    if all(isnan.(initial_conditions.initial_viscosity))
        default_viscosity = params.default_viscosity
    #  @info "Did not find a specified initial viscosity, reverting to default value specified in params ($default_viscosity Pa s everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_viscosity =  default_viscosity*ones(grid.nx, grid.ny, grid.nσ)
    else
        initial_viscosity = initial_conditions.initial_viscosity
     end

    if all(isnan.(initial_conditions.initial_temperature))
        default_temperature = params.default_temperature
        #@info "Did not find a specified initial temperature, reverting to default value specified in params ($default_temperature K everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_temperature =  default_temperature*ones(grid.nx, grid.ny, grid.nσ)
    else
        initial_temperature = initial_conditions.initial_temperature
    end

    if all(isnan.(initial_conditions.initial_damage))
        default_damage = params.default_damage
    #  @info "Did not find a specified initial damage field, reverting to default value specified in params ($default_damage everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_damage =  default_damage*ones(grid.nx, grid.ny, grid.nσ)
    else
        initial_damage = initial_conditions.initial_damage
    end

    if all(isnan.(initial_conditions.initial_strain_history))
        default_strain_history = params.default_strain_history
        #@info "Did not find a specified initial tensile strain history  field, reverting to default value specified in params ($default_strain_history everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_strain_history =  default_strain_history*ones(grid.nx, grid.ny, grid.nσ)
    else
        initial_strain_history = initial_conditions.initial_strain_history
    end

    if all(isnan.(initial_conditions.initial_basal_water_thickness))
        default_basal_water_thickness = params.basal_water_thickness
        #@info "Did not find a specified initial basal water thickness field, reverting to default value specified in params ($default_basal_water_thickness everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_basal_water_thickness =  default_basal_water_thickness*ones(grid.nx, grid.ny)
    else
        initial_basal_water_thickness = initial_conditions.initial_basal_water_thickness
    end

    if all(isnan.(initial_conditions.initial_hydraulic_potential_b))
        hydraulic_potential_b = params.hydraulic_potential_b
        #@info "Did not find a specified initial hydraulic_potential field, reverting to default value specified in params ($hydraulic_potential_b everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_hydraulic_potential_b =  hydraulic_potential_b.*ones(grid.nx, grid.ny)
    else
        initial_hydraulic_potential_b = initial_conditions.initial_hydraulic_potential_b
    end

    if all(isnan.(initial_conditions.initial_effective_pressure))
        default_effective_pressure = params.effective_pressure
        #@info "Did not find a specified initial effective pressure field, reverting to default value specified in params ($effective_pressure everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_effective_pressure =  default_effective_pressure.*ones(grid.nx, grid.ny)
    else
        initial_effective_pressure = initial_conditions.initial_effective_pressure
    end

    if all(isnan.(initial_conditions.initial_θ_ave))
        default_θ_ave = params.default_temperature_ave
        #@info "Did not find a specified initial depth-averaged temperature field, reverting to default value specified in params ($default_θ_ave everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_θ_ave =  default_θ_ave*ones(grid.nx, grid.ny)
    else
        initial_θ_ave = initial_conditions.initial_θ_ave
    end

    if all(isnan.(initial_conditions.initial_preBfactor))
        default_preBfactor = params.default_preBfactor
        #@info "Did not find a specified initial depth-averaged temperature field, reverting to default value specified in params ($default_θ_ave everywhere)...\n...If you have set niter0 > 0 without invoking the update flag, you can ignore this message"
        initial_preBfactor =  default_preBfactor*ones(grid.nx, grid.ny)
    else
        initial_preBfactor = initial_conditions.initial_preBfactor
    end

    #check sizes are compatible
    (size(initial_thickness) == (grid.nx, grid.ny)) || throw(DimensionMismatch("Initial thickness field is not compatible with grid size. Input thickess field is has size $(size(initial_conditions.initial_thickness)), which must match horizontal grid size ($(grid.nx) x $(grid.ny))"))
    (size(initial_grounded_fraction) == (grid.nx, grid.ny)) || throw(DimensionMismatch("Initial grounded fraction field is not compatible with grid size. Input grounded fraction field is has size $(size(initial_conditions.initial_grounded_fraction)), which must match horizontal grid size ($(grid.nx) x $(grid.ny))"))
    (size(initial_u_veloc) == (grid.nx+1, grid.ny)) || throw(DimensionMismatch("Initial u_velocity field is not compatible with grid size. Input u_velocity field has size $(size(initial_conditions.initial_u_veloc)), which must match horizontal grid size ($(grid.nx+1) x $(grid.ny))"))
    (size(initial_v_veloc) == (grid.nx, grid.ny+1)) || throw(DimensionMismatch("Initial v_velocity field is not compatible with grid size. Input v_velocity field has size $(size(initial_conditions.initial_v_veloc)), which must match horizontal grid size ($(grid.nx) x $(grid.ny+1))"))
    (size(initial_temperature) == (grid.nx, grid.ny, grid.nσ)) || throw(DimensionMismatch("Initial temperature field is not compatible with grid size. Input temperature field is has size $(size(initial_conditions.initial_temperature)), which must match 3D grid size ($(grid.nx), $(grid.ny), $(grid.nσ))"))
    (size(initial_viscosity) == (grid.nx, grid.ny, grid.nσ)) || throw(DimensionMismatch("Initial viscosity field is not compatible with grid size. Input viscosity field is has size $(size(initial_conditions.initial_viscosity)), which must match 3D grid size ($(grid.nx), $(grid.ny), $(grid.nσ))"))
    (size(initial_damage) == (grid.nx, grid.ny, grid.nσ)) || throw(DimensionMismatch("Initial damage field is not compatible with grid size.Input damage field is has size $(size(initial_conditions.initial_damage)), which must match 3D grid size ($(grid.nx), $(grid.ny), $(grid.nσ))"))
    (size(initial_strain_history) == (grid.nx, grid.ny, grid.nσ)) || throw(DimensionMismatch("Initial tensile strain history field is not compatible with grid size.Input tensile strain history field has size $(size(initial_conditions.initial_strain_history)), which must match 3D grid size ($(grid.nx), $(grid.ny), $(grid.nσ))"))
    (size(initial_basal_water_thickness) == (grid.nx, grid.ny)) || throw(DimensionMismatch("Initial basal water thickness field is not compatible with grid size. Input basal water thickness field has size $(size(initial_conditions.initial_basal_water_thickness)), which must match horizontal grid size ($(grid.nx) x $(grid.ny))"))
    (size(initial_hydraulic_potential_b) == (grid.nx, grid.ny)) || throw(DimensionMismatch("Initial hydraulic_potential field is not compatible with grid size. Input hydraulic_potential field has size $(size(initial_conditions.initial_hydraulic_potential_b)), which must match horizontal grid size ($(grid.nx) x $(grid.ny))"))
    (size(initial_effective_pressure) == (grid.nx, grid.ny)) || throw(DimensionMismatch("Initial effective pressure field is not compatible with grid size. Input effective pressure field has size $(size(initial_conditions.initial_effective_pressure)), which must match horizontal grid size ($(grid.nx) x $(grid.ny))"))
    (size(initial_θ_ave) == (grid.nx, grid.ny)) || throw(DimensionMismatch("Initial depth-averaged temperature field is not compatible with grid size. Input depth-averaged temperature field has size $(size(initial_conditions.initial_θ_ave)), which must match horizontal grid size ($(grid.nx) x $(grid.ny))"))
    (size(initial_preBfactor) == (grid.nx, grid.ny)) || throw(DimensionMismatch("Initial preBfactor is not compatible with grid size. Input preBfactor field has size $(size(initial_conditions.initial_preBfactor)), which must match horizontal grid size ($(grid.nx) x $(grid.ny))"))

    return     InitialConditions(
        initial_thickness,
        initial_grounded_fraction,
        initial_u_veloc,
        initial_v_veloc,
        initial_viscosity,
        initial_temperature,
        initial_damage,
        initial_strain_history,
        initial_basal_water_thickness,
        initial_hydraulic_potential_b,
        initial_effective_pressure,
        initial_θ_ave,
        initial_preBfactor
    )
end    

function reconstruct_on_subdomain(initial_conditions::InitialConditions, grid::Grid, subdomain::NTuple{4,<: Integer}) 
    
    x_start,x_end,y_start,y_end = subdomain
    u_grid_size, v_grid_size = (grid.nx+1, grid.ny), (grid.nx, grid.ny+1)
    
    return InitialConditions(
        initial_thickness = 
        size(initial_conditions.initial_thickness) == size(grid)[1:2] ? 
            initial_conditions.initial_thickness[x_start:x_end, y_start:y_end] : 
            initial_conditions.initial_thickness,
        initial_grounded_fraction = 
        size(initial_conditions.initial_grounded_fraction) == size(grid)[1:2] ? 
            initial_conditions.initial_grounded_fraction[x_start:x_end, y_start:y_end] : 
            initial_conditions.initial_grounded_fraction,
        initial_u_veloc = 
        size(initial_conditions.initial_u_veloc) == u_grid_size ? 
            initial_conditions.initial_u_veloc[x_start:x_end+1, y_start:y_end] : 
            initial_conditions.initial_u_veloc,
        initial_v_veloc = 
        size(initial_conditions.initial_v_veloc) == v_grid_size ? 
            initial_conditions.initial_v_veloc[x_start:x_end, y_start:y_end+1] : 
            initial_conditions.initial_v_veloc,
        initial_viscosity = 
        size(initial_conditions.initial_viscosity) == size(grid) ? 
            initial_conditions.initial_viscosity[x_start:x_end, y_start:y_end, :] : 
            initial_conditions.initial_viscosity,
        initial_temperature = 
        size(initial_conditions.initial_temperature) == size(grid) ? 
            initial_conditions.initial_temperature[x_start:x_end, y_start:y_end, :] :
            initial_conditions.initial_temperature,
        initial_damage = 
        size(initial_conditions.initial_damage) == size(grid) ? 
            initial_conditions.initial_damage[x_start:x_end, y_start:y_end, :] : 
            initial_conditions.initial_damage,
        initial_strain_history =
        size(initial_conditions.initial_strain_history) == size(grid) ? 
            initial_conditions.initial_strain_history[x_start:x_end, y_start:y_end, :] : 
            initial_conditions.initial_strain_history,
        initial_basal_water_thickness =
        size(initial_conditions.initial_basal_water_thickness) == size(grid)[1:2] ? 
            initial_conditions.initial_basal_water_thickness[x_start:x_end, y_start:y_end] : 
            initial_conditions.initial_basal_water_thickness,
        initial_hydraulic_potential_b =
        size(initial_conditions.initial_hydraulic_potential_b) == size(grid)[1:2] ? 
            initial_conditions.initial_hydraulic_potential_b[x_start:x_end, y_start:y_end] : 
            initial_conditions.initial_hydraulic_potential_b,    
        initial_effective_pressure =
        size(initial_conditions.initial_effective_pressure) == size(grid)[1:2] ? 
            initial_conditions.initial_effective_pressure[x_start:x_end, y_start:y_end] : 
            initial_conditions.initial_effective_pressure,
        initial_θ_ave =
        size(initial_conditions.initial_θ_ave) == size(grid)[1:2] ? 
            initial_conditions.initial_θ_ave[x_start:x_end, y_start:y_end] : 
            initial_conditions.initial_θ_ave,
        initial_preBfactor =
        size(initial_conditions.initial_preBfactor) == size(grid)[1:2] ? 
            initial_conditions.initial_preBfactor[x_start:x_end, y_start:y_end] : 
            initial_conditions.initial_preBfactor
    )

end

# Maps each InitialConditions field to a path on GridField.
# MPI gather/scatter uses the same list so it won't diverge from a serial run.
const INITIAL_CONDITION_FIELD_PATHS = (
    :initial_thickness => (:gh, :h),
    :initial_grounded_fraction => (:gh, :grounded_fraction),
    :initial_u_veloc => (:gu, :u),
    :initial_v_veloc => (:gv, :v),
    :initial_viscosity => (:g3d, :η),
    :initial_temperature => (:g3d, :θ),
    :initial_damage => (:g3d, :Φ),
    :initial_strain_history => (:g3d, :strain_history),
    :initial_basal_water_thickness => (:gh, :basal_water_thickness),
    :initial_hydraulic_potential_b => (:gh, :hydraulic_potential_b),
    :initial_effective_pressure => (:gh, :effective_pressure),
    :initial_θ_ave => (:gh, :θ_ave),
    :initial_preBfactor => (:gh, :preBfactor),
)

const CHECKPOINT_EXTRA_FIELD_PATHS = ((:gh, :β),)
const CHECKPOINT_FIELD_PATHS = (map(last, INITIAL_CONDITION_FIELD_PATHS)..., CHECKPOINT_EXTRA_FIELD_PATHS...)

function gridfield_get(fields, path::Tuple)
    value = fields
    for name in path
        value = getproperty(value, name)
    end
    return value
end

"""Copy restart arrays from a `GridField` into an `InitialConditions` object."""
function initial_conditions_from_fields(fields)
    kwargs = Pair{Symbol, Any}[name => copy(gridfield_get(fields, path)) for (name, path) in INITIAL_CONDITION_FIELD_PATHS]
    return InitialConditions(; kwargs...)
end

"""
Copy restart arrays from `initial_conditions` into a live `GridField`.

`names` selects which IC fields to copy. A size mismatch is an error, not a skip.
"""
function copy_initial_conditions_to_fields!(
    fields,
    initial_conditions,
    names = fieldnames(typeof(initial_conditions)),
)
    for (name, path) in INITIAL_CONDITION_FIELD_PATHS
        name in names || continue
        src = getfield(initial_conditions, name)
        dst = gridfield_get(fields, path)
        if size(src) != size(dst)
            error("Checkpoint $name has size $(size(src)), live field has $(size(dst))")
        end
        dst .= src
    end
    return nothing
end

