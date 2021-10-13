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
Includes an option for not evolving shelf thickness
"""
function update_thickness!(simulation)
@unpack model,timestepping_params=simulation
hUpdate=zeros(model.grid.nx,model.grid.ny)
aground=zeros(model.grid.nx,model.grid.ny)
hUpdate[model.fields.gh.mask]=max.(model.params.minimum_thickness .- model.fields.gh.h[model.fields.gh.mask],timestepping_params.dt*model.fields.gh.dhdt[model.fields.gh.mask])
#Specify whether to evolve the shelves:
if !model.params.evolveShelves
    hUpdate[model.fields.gh.mask]=max.(model.params.smallHAF.-(model.params.density_ocean./model.params.density_ice).*model.fields.gh.b[model.fields.gh.mask].-model.fields.gh.h[model.fields.gh.mask],hUpdate[model.fields.gh.mask])
    aground=(model.fields.gh.haf.>=0)
    wc=[1 1 1; 1 1 1; 1 1 1]
    w=centered(wc)
    nearfloat_mask=imfilter(model.fields.gh.mask.&.!aground,reflect(w),Fill(0,w))
    nearfloat_mask=iszero.(iszero.(nearfloat_mask))
    hUpdate[nearfloat_mask].=0
  end
model.fields.gh.h[model.fields.gh.mask] = model.fields.gh.h[model.fields.gh.mask] .+ hUpdate[model.fields.gh.mask]
return simulation
end

"""
    update_clock!(simulation::AbstractSimulation)

Update the simulation clock
"""
function update_clock!(simulation)
    @unpack clock,timestepping_params=simulation
    clock.n_iter += 1
    clock.time += timestepping_params.dt
    println("The time is ", clock.time)
    return simulation
end


"""
    run_simulation(simulation)
Perform the simulation specified by the simulation
"""
function run_simulation!(simulation)
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
            println(simulation.clock.n_iter)
            println("outputting at timestep number $(simulation.clock.n_iter)")

        end

        #check the dump velocity flag at the final timestep
        if (i == timestepping_params.n_iter_total) && output_params.dump_vel
            write_vel(simulation)
        end
    end

    #zip the simulation output (no zipping handled by zip_output)
    zip_output(simulation)
        
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
    
    u_out=model.fields.gu.u[1:end-1,:]
    v_out=model.fields.gv.v[:,1:end-1]

    u_out .= hton.(u_out)
    v_out .= hton.(v_out)

    ufileID =  open(uVel_file_string,"w")
      write(ufileID, u_out[:,:])
    close(ufileID) 
    vfileID =  open(vVel_file_string,"w")
    write(vfileID, v_out[:,:])
    close(vfileID)   

 end 