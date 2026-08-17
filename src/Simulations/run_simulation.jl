export run_simulation!, timestep!, update_clock!, update_thickness!, update_model_climate_forcing!

using WAVI
import WAVI: AbstractModel, AbstractSimulation
using WAVI.Outputs: write_outputs, zip_output
using WAVI.Processes: update_state!, update_state_novelocity!

using WAVI.MeltRates
using WAVI.SurfaceMassBalance
using WAVI.Fracture
using WAVI.SlidingLaw
using WAVI.BasalHydrology
using WAVI.ThermoDynamics


"""
    timestep!(model, output_params, clock, timestepping_params)

Perform one timestep of the simulation. 
"""
function timestep!(model::AbstractModel{T,N,S},
                   timestepping_params::TimesteppingParams,
                   output_params::OutputParams,
                   clock::Clock) where {T,N,S}

 
    if mod(clock.n_iter, timestepping_params.ntimesteps_climate_forcing_update) == 0 #if we're at a number of timesteps to update the forcing files, do that
        update_model_climate_forcing!(model, clock)
    end

    if mod(clock.n_iter, timestepping_params.ntimesteps_velocity_update) == 0 #if we're on a velocity + thickness update output step
        update_state!(model, clock)
        if timestepping_params.verbose
            println("Completed main timestep update), t = ", clock.time)
        end
    else #just a thickness update step
        update_state_novelocity!(model, clock)
        if timestepping_params.verbose
            println("Completed sub-timestep update, t = ", clock.time)
        end
    end

    #write solution if at the first timestep (hack for https://github.com/RJArthern/WAVI.jl/issues/46 until synchronicity is fixed)
    # Have made the interface consistent
    # Have also removed the dependence on individual calls
    if (output_params.output_start) && (clock.n_iter == 0)
        write_outputs(model, timestepping_params, output_params, clock)
    end


    if timestepping_params.step_thickness
        update_thickness!(model, timestepping_params)
    end

    if timestepping_params.verbose
        println("Completed thickness update), t = ", clock.time)
    end

    update_clock!(clock, timestepping_params)

    write_outputs(model, timestepping_params, output_params, clock)

end
timestep!(s::AbstractSimulation) = timestep!(s.model, s.timestepping_params, s.output_params, s.clock)

"""
update_thickness!(model::AbstractModel)

Update thickness using rate of change of thickness and apply minimum thickness constraint. Includes an option for not evolving shelf thickness.
"""
function update_thickness!(model::AbstractModel, timestepping_params)
    hUpdate = zeros(model.grid.nx,model.grid.ny)
    aground = zeros(model.grid.nx,model.grid.ny)
    hUpdate[model.fields.gh.mask] = max.(
        model.params.minimum_thickness .- model.fields.gh.h[model.fields.gh.mask],
        timestepping_params.dt * model.fields.gh.dhdt[model.fields.gh.mask])
    
    #Specify whether to evolve the shelves:
    if !model.params.evolveShelves
        hUpdate[model.fields.gh.mask] = max.(
            model.params.smallHAF .- (
                model.params.density_ocean ./ model.params.density_ice
            ) .* model.fields.gh.b[model.fields.gh.mask] .- model.fields.gh.h[model.fields.gh.mask], 
            hUpdate[model.fields.gh.mask])
        aground=(model.fields.gh.haf.>=0)
        wc=[1 1 1; 1 1 1; 1 1 1]
        w=centered(wc)
        nearfloat_mask = imfilter(model.fields.gh.mask.&.!aground,reflect(w),Fill(0,w))
        nearfloat_mask = iszero.(iszero.(nearfloat_mask))
        hUpdate[nearfloat_mask].=0
    end
    hUpdate[model.fields.gh.h_isfixed] .= 0
    model.fields.gh.h[model.fields.gh.mask] = model.fields.gh.h[model.fields.gh.mask] .+ hUpdate[model.fields.gh.mask]
end
update_thickness!(simulation::Simulation) = update_thickness!(s.model, s.timestepping_params)


"""
    update_clock!(simulation::AbstractSimulation)

Update the simulation clock
"""
function update_clock!(clock::Clock, timestepping_params::TimesteppingParams)
    clock.n_iter += 1
    clock.time += timestepping_params.dt
end
update_clock!(s::Simulation) = update_clock!(s.clock, s.timestepping_params)

"""
    run_simulation!(simulation)
    
Perform the simulation specified by the simulation
"""
function run_simulation!(model::AbstractModel{T,N,S}, 
                         timestepping_params::TimesteppingParams, 
                         output_params::OutputParams,
                         clock::Clock) where {T,N,S}
    spec = model.spec
    spec_label = hasproperty(spec, :ngridsx) ?
        "$(nameof(typeof(spec))) $(spec.ngridsx)x$(spec.ngridsy)" :
        string(nameof(typeof(spec)))
    @info "$(spec_label) grid $(model.grid.nx)x$(model.grid.ny)"
    for i = (clock.n_iter+1):timestepping_params.n_iter_total
        @info "Running iteration $(clock.n_iter)/$(timestepping_params.n_iter_total)"
        if clock.n_iter == 0
            @info "First velocity solve compiles kernels and may take a while"
        end
        timestep!(model, timestepping_params, output_params, clock)
    end

    #zip the simulation output (no zipping handled by zip_output)
    zip_output(model, output_params)
end
run_simulation!(s::Simulation) = run_simulation!(s.model, s.timestepping_params, s.output_params, s.clock)


function update_model_climate_forcing!(model::AbstractModel, clock)
        @unpack shelf_melt_rate, surface_mass_balance, fracture, sliding_law, basal_hydrology, thermo_dynamics, grid = model

        update_climate_forcing!(shelf_melt_rate, grid, clock)
        update_climate_forcing!(surface_mass_balance, grid, clock)
        update_climate_forcing!(fracture, grid, clock)
        update_climate_forcing!(sliding_law, grid, clock)
        update_climate_forcing!(basal_hydrology, grid, clock)
        update_climate_forcing!(thermo_dynamics, grid, clock)

    return nothing
end