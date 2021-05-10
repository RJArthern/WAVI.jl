"""
    timestep!(simulation)

Perform one timestep of the simulation
"""
function timestep!(simulation)
    @unpack model,timestepping_params = simulation
    update_state!(model)
    update_thickness!(simulation)
    update_clock!(simulation)
end

"""
update_thickness!(model::AbstractModel)

Update thickness using rate of change of thickness and apply minimum thickness constraint.
"""
function update_thickness!(simulation::AbstractSimulation)
@unpack model,timestepping_params=simulation
onesvec=ones(model.grid.nx*model.grid.ny)
model.gh.h[model.gh.mask] = model.gh.h[model.gh.mask] .+ max.(model.params.minimum_thickness .- model.gh.h[model.gh.mask],timestepping_params.dt*model.gh.dhdt[model.gh.mask])
println("hello")
return simulation
end

"""
    update_clock!(simulation::AbstractSimulation)

Update the simulation clock
"""
function update_clock!(simulation::AbstractSimulation)
    @unpack clock,timestepping_params=simulation
    clock.n_iter += 1
    clock.time += timestepping_params.dt
    return simulation
end


"""
    run_simulation(simulation)
Perform the simulation specified by the simulation
"""
function run_simulation(simulation::Simulation)
    @unpack model, timestepping_params, output_params = simulation

    ~(timestepping_params === nothing) || error("Must pass timestepping params")

    if timestepping_params.n_iter0 == 0 #start a fresh run
        #check that grid and bed have been passed
        ~(grid === nothing) || error("Must pass a grid if starting a fresh run")
        ~(bed_elevation === nothing) || error("Must pass a bed (array or function) if starting a fresh run")

        println("Starting clean wavi simulation")


        #if no parameters have been passed, construct defaults
        if (params === nothing); params = Params(); end 
        if (solver_params === nothing); solver_params = SolverParams(); end 
        if (initial_conditions === nothing); initial_conditions = InitialConditions(); end #don't worry about strange defaults here, these will be picked up by start
        if (timestepping_params === nothing); timestepping_params = TimesteppingParams(); end 
        if (output == nothing); output = Output(); end
        wavi = start(grid = grid, 
                    bed_elevation = bed_elevation,
                    params = params, 
                    solver_params = solver_params,
                    initial_conditions = initial_conditions,
                    timestepping_params = timestepping_params,
                    output = output)

        #initialize things
        chkpt_tag = "A" #initialize the checkpoint tag
        println("running simulation...")

        #timestepping loop
        for i = 1:wavi.timestepping_params.n_iter_total
            run!(wavi)

            #check if we have hit a temporary checkpoint
            if mod(i,wavi.timestepping_params.n_iter_chkpt) == 0
                #output a temporary checkpoint
                fname = string("Chkpt",chkpt_tag, ".jld2")
                @save fname wavi
                chkpt_tag = (chkpt_tag == "A") ? "B" : "A"
                println("making temporary checkpoint at iteration number $(wavi.clock.n_iter)")
            end

            #check if we have hit a permanent checkpoint
            if mod(i,wavi.timestepping_params.n_iter_pchkpt) == 0
                #output a permanent checkpoint
                n_iter_string =  lpad(wavi.clock.n_iter, 10, "0"); #filename as a string with 10 digits
                fname = string("PChkpt_",n_iter_string, ".jld2")
                @save fname wavi
                println("making permanent checkpoint at iteration number $(wavi.clock.n_iter)")
            end

            #check if we have hit an output timestep
            if mod(i,wavi.output.n_iter_out) == 0
                write_output(wavi)
            end
        end
        

    else #look for a pickup
        n_iter_string =  lpad(timestepping_params.n_iter0, 10, "0"); #filename as a string with 10 digits
        try 
            @load string("PChkpt_",n_iter_string, ".jld2") wavi
            println("Pickup successful")
        catch 
            println("Pickup error, terminating run")
        end

        #update the parameters of those that have been specified the flag is specified

        #continue with the run
        chkpt_tag = "A" #initialize the checkpoint tag
        for i = (timestepping_params.n_iter0+1):timestepping_params.n_iter_total
            run!(wavi)
            if mod(i,timestepping_params.n_iter_chkpt) == 0
                #output a temporary checkpoint
                println("making temporary checkpoint at iteration number $(wavi.clock.n_iter)")
                fname = string("Chkpt",chkpt_tag, ".jld2")
                @save fname wavi
                chkpt_tag = (chkpt_tag == "A") ? "B" : "A"
            end
            if mod(i,timestepping_params.n_iter_pchkpt) ==0 
                #output a permanent checkpoint
                println("making permanent checkpoint at iteration number $(wavi.clock.n_iter)")
                n_iter_string =  lpad(wavi.clock.n_iter, 10, "0"); #filename as a string with 10 digits
                fname = string("PChkpt_",n_iter_string, ".jld2")
                @save fname wavi

            end
        end
    end
    return wavi
end

