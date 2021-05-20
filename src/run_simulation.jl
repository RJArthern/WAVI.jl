"""
    timestep!(simulation)

Perform one timestep of the simulation
"""
function timestep!(simulation)
    @unpack model,timestepping_params = simulation
    update_state!(model)
    if timestepping_params.step_thickness
        update_thickness!(simulation)
    end
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
function run_simulation!(simulation::Simulation)
    @unpack model, timestepping_params, output_params = simulation
    chkpt_tag = "A"
    for i = (simulation.clock.n_iter+1):timestepping_params.n_iter_total
        timestep!(simulation)

        #check if we have hit a temporary checkpoint
        if mod(i,timestepping_params.n_iter_chkpt) == 0
            #output a temporary checkpoint
            fname = string("Chkpt",chkpt_tag, ".jld2")
            @save fname simulation
            chkpt_tag = (chkpt_tag == "A") ? "B" : "A"
            println("making temporary checkpoint at timestep number $(simulation.clock.n_iter)")
        end

        #check if we have hit a permanent checkpoint
        if mod(i,simulation.timestepping_params.n_iter_pchkpt) == 0
            #output a permanent checkpoint
            n_iter_string =  lpad(simulation.clock.n_iter, 10, "0"); #filename as a string with 10 digits
            fname = string("PChkpt_",n_iter_string, ".jld2")
            @save fname simulation
            println("making permanent checkpoint at timestep number $(simulation.clock.n_iter)")
        end

        #check if we have hit an output timestep
        if mod(i,simulation.output_params.n_iter_out) == 0
            write_output(simulation)
            println("outputting at timestep number $(simulation.clock.n_iter)")

        end

        #check the dump velocity flag at the final timestep
        if (i == timestepping_params.n_iter_total) && output_params.dump_vel
            write_vel(simulation)
        end
    end
        
    return simulation
end

"""
    function write_vel(simulation)

Write the velocity at the the final timestep of the simulation (used in the coupled wavi-mitgcm model to communicate with streamice)
"""
function write_vel(simulation::Simulation)
    @unpack model = simulation  
    uVel_file_string = string(simulation.output_params.prefix,  "_U.bin")
    vVel_file_string = string(simulation.output_params.prefix,  "_V.bin")
    
    u_out=model.gh.u 
    v_out=model.gh.v

    u_out .= hton.(u_out)
    v_out .= hton.(v_out)

    ufileID =  open(uVel_file_string,"w")
      write(ufileID, u_out[:,:])
    close(ufileID) 
    vfileID =  open(vVel_file_string,"w")
    write(vfileID, v_out[:,:])
    close(vfileID)   

 end 