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
function update_thickness!(simulation)
@unpack model,timestepping_params=simulation
onesvec=ones(model.grid.nx*model.grid.ny)
    
if model.grid.Cxl > 1 || model.grid.Cxu <  model.grid.nx || model.grid.Cyl > 1 || model.grid.Cyu < model.grid.ny
 
 thick_nochild=model.fields.gh.h[model.grid.Cxl:model.grid.Cxu,model.grid.Cyl:model.grid.Cyu]  

end        
    
model.fields.gh.h[model.fields.gh.mask] = model.fields.gh.h[model.fields.gh.mask] .+ max.(model.params.minimum_thickness .- model.fields.gh.h[model.fields.gh.mask],timestepping_params.dt*model.fields.gh.dhdt[model.fields.gh.mask])

if model.grid.Cxl > 1 || model.grid.Cxu <  model.grid.nx || model.grid.Cyl > 1 || model.grid.Cyu < model.grid.ny
model.fields.gh.h[model.grid.Cxl:model.grid.Cxu,model.grid.Cyl:model.grid.Cyu]=thick_nochild[:,:]
end

        
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
    return simulation
end


"""
    run_simulation(simulation)
Perform the simulation specified by the simulation
"""
function run_simulation!(simulation)
    @unpack model, timestepping_params, output_params = simulation
    chkpt_tag = "A"
    if output_params.dump_vel
       h_out_line_w = zeros(model.grid.Cyu - model.grid.Cyl +1) 
       h_out_line_e = zeros(model.grid.Cyu - model.grid.Cyl +1)
       h_out_line_s = zeros(model.grid.Cxu - model.grid.Cxl +1) 
       h_out_line_n = zeros(model.grid.Cxu - model.grid.Cxl +1) 
    end
    for i = (simulation.clock.n_iter+1):timestepping_params.n_iter_total
        timestep!(simulation)
        #if output_params.dump_vel
        # if output_params.PC_west
        #  if model.grid.Cxl > 1
        #   h_out_line_w = h_out_line_w + model.fields.gh.h[model.grid.Cxl-1,model.grid.Cyl:model.grid.Cyu]
        #    if (i == timestepping_params.n_iter_total)
        #     h_out_line_w= h_out_line_w ./ (timestepping_params.n_iter_total- timestepping_params.niter0)
        #    end
        #  end
        # end
        # if output_params.PC_east
        #  if model.grid.Cxu < model.grid.nx
        #  h_out_line_e = h_out_line_e + model.fields.gh.h[model.grid.Cxu + 1,model.grid.Cyl:model.grid.Cyu]
        #   if (i == timestepping_params.n_iter_total)
        #    h_out_line_e= h_out_line_e ./ (timestepping_params.n_iter_total- timestepping_params.niter0)
        #   end
        #  end
        # end
        # if output_params.PC_south
        #  if model.grid.Cyl > 1
        #  h_out_line_s = h_out_line_s + model.fields.gh.h[model.grid.Cxl:model.grid.Cxu,model.grid.Cyl-1]
        #   if (i == timestepping_params.n_iter_total)
        #   h_out_line_s= h_out_line_s ./ (timestepping_params.n_iter_total- timestepping_params.niter0)
        #   end
        #  end
        # end
        # if output_params.PC_north
        #  if model.grid.Cyu < model.grid.ny
        #  h_out_line_n = h_out_line_n + model.fields.gh.h[model.grid.Cxl:model.grid.Cxu,model.grid.Cyu+1]
        #   if (i == timestepping_params.n_iter_total)
        #   h_out_line_n= h_out_line_n ./ (timestepping_params.n_iter_total- timestepping_params.niter0)
        #   end
        #  end
        # end
        #end
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
         if output_params.PC_west
          if model.grid.Cxl > 1
           h_out_line_w = model.fields.gh.h[model.grid.Cxl-1,model.grid.Cyl:model.grid.Cyu]
          end
         end
         if output_params.PC_east
          if model.grid.Cxu < model.grid.nx
          h_out_line_e = model.fields.gh.h[model.grid.Cxu + 1,model.grid.Cyl:model.grid.Cyu]
          end
         end  
         if output_params.PC_south
          if model.grid.Cyl > 1
          h_out_line_s = model.fields.gh.h[model.grid.Cxl:model.grid.Cxu,model.grid.Cyl-1]
          end
         end
         if output_params.PC_north
          if model.grid.Cyu < model.grid.ny
          h_out_line_n = model.fields.gh.h[model.grid.Cxl:model.grid.Cxu,model.grid.Cyu+1]
          end
         end
            write_vel(simulation,h_out_line_w,h_out_line_e,h_out_line_n,h_out_line_s)
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
function write_vel(simulation::Simulation,h_out_line_w,h_out_line_e,h_out_line_n,h_out_line_s)
    @unpack model, output_params = simulation  
    x_w=0  
    x_e=0     
    y_s=0
    y_n=0
    if model.grid.Cxl > 1 || model.grid.Cxu <  model.grid.nx || model.grid.Cyl > 1 || model.grid.Cyu < model.grid.ny  
     if output_params.PC_west
      if model.grid.Cxl > 1 
       x_w=x_w+1
      end
     end
     if output_params.PC_east
      if model.grid.Cxu <  model.grid.nx
       x_e=x_e+1
      end
     end
     if output_params.PC_south
      if model.grid.Cyl > 1 
       y_s=y_s+1
      end
     end
     if output_params.PC_north
      if model.grid.Cyu <  model.grid.ny
       y_n=y_n+1
      end
     end
    end
    
    uVel_file_string = string(simulation.output_params.prefix,  "_U.bin")
    vVel_file_string = string(simulation.output_params.prefix,  "_V.bin")
    
    u_out = zeros(model.grid.Cxu - model.grid.Cxl + 1 + 2 + x_e +x_w,model.grid.Cyu - model.grid.Cyl +1 +2 + y_s +y_n)
    v_out = zeros(model.grid.Cxu - model.grid.Cxl + 1 + 2 + x_e +x_w,model.grid.Cyu - model.grid.Cyl +1 +2 + y_s +y_n)
    
    u_out[2:end-1,2:end-1]=model.fields.gu.u[model.grid.Cxl-x_w:model.grid.Cxu+x_e,model.grid.Cyl-y_s:model.grid.Cyu+y_n]
    v_out[2:end-1,2:end-1]=model.fields.gv.v[model.grid.Cxl-x_w:model.grid.Cxu+x_e,model.grid.Cyl-y_s:model.grid.Cyu+y_n]

    u_out .= hton.(u_out)
    v_out .= hton.(v_out)

    ufileID =  open(uVel_file_string,"w")
      write(ufileID, u_out[:,:])
    close(ufileID) 
    vfileID =  open(vVel_file_string,"w")
    write(vfileID, v_out[:,:])
    close(vfileID)   
    
    if model.grid.Cxl > 1 || model.grid.Cxu <  model.grid.nx || model.grid.Cyl > 1 || model.grid.Cyu < model.grid.ny 
            
            
     h_out_b = zeros(model.grid.Cxu - model.grid.Cxl + 1 + 2 + x_e +x_w ,model.grid.Cyu - model.grid.Cyl +1 + 2 + y_s + y_n)
     if output_params.PC_west       
      if model.grid.Cxl > 1
                h_out_b[2,y_s+2:end-y_n-1] .= h_out_line_w[:]
      end
     end
     if output_params.PC_east
      if model.grid.Cxu < model.grid.nx
                h_out_b[end-1,y_s+2:end-y_n-1] .= h_out_line_e[:]
      end
     end
     if output_params.PC_south
      if model.grid.Cyl > 1
                h_out_b[x_w+2:end-x_e-1,2] .= h_out_line_s[:]
      end
     end
     if output_params.PC_north
      if model.grid.Cyu < model.grid.ny
                h_out_b[x_w+2:end-x_e-1,end-1] .= h_out_line_n[:]
      end
     end
     h_out_b .= hton.(h_out_b)
     
     clock_time_s=Int(simulation.clock.time*(3600*24*365))
        
     hb_file_string = string(simulation.output_params.prefix, "_Hb", lpad(clock_time_s, 10,"0") ,".bin")
    
     hbfileID =  open(hb_file_string,"w")
     write(hbfileID, h_out_b[:,:])
     close(hbfileID) 
        
    end
 end 
