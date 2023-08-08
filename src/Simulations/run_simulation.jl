"""
    timestep!(simulation)

Perform one timestep of the simulation
"""
function timestep!(simulation)
    @unpack model,timestepping_params, output_params, clock = simulation
    update_state!(model, clock)
    #write solution if at the first timestep (hack for https://github.com/RJArthern/WAVI.jl/issues/46 until synchronicity is fixed)
    if (output_params.output_start) && (simulation.clock.n_iter == 0)
        write_output(simulation)
    end
    if timestepping_params.step_thickness
        update_thickness!(simulation)
    end
    update_clock!(simulation)
end

"""
update_thickness!(model::AbstractModel)

Update thickness using rate of change of thickness and apply minimum thickness constraint. Includes an option for not evolving shelf thickness.
"""
function update_thickness!(simulation)
@unpack model,timestepping_params=simulation
onesvec=ones(model.grid.nx*model.grid.ny)
hUpdate=zeros(model.grid.nx,model.grid.ny)
aground=zeros(model.grid.nx,model.grid.ny)
    
if model.grid.Cxl > 1 || model.grid.Cxu <  model.grid.nx || model.grid.Cyl > 1 || model.grid.Cyu < model.grid.ny
 
 thick_nochild=model.fields.gh.h[model.grid.Cxl:model.grid.Cxu,model.grid.Cyl:model.grid.Cyu]  

end        
        
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
hUpdate[model.fields.gh.h_isfixed] .= 0
model.fields.gh.h[model.fields.gh.mask] = model.fields.gh.h[model.fields.gh.mask] .+ hUpdate[model.fields.gh.mask]

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
       h_out_line_w = zeros(model.grid.Cyu - model.grid.Cyl +3) 
       h_out_line_e = zeros(model.grid.Cyu - model.grid.Cyl +3)
       h_out_line_s = zeros(model.grid.Cxu - model.grid.Cxl +3) 
       h_out_line_n = zeros(model.grid.Cxu - model.grid.Cxl +3) 
    end
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
        end

        #check the dump velocity flag at the final timestep
        if (i == timestepping_params.n_iter_total) && output_params.dump_vel
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
            
         if output_params.PC_west
          if model.grid.Cxl > 1
           h_out_line_w = model.fields.gh.h[model.grid.Cxl-1,model.grid.Cyl-y_s:model.grid.Cyu+y_n]
          end
         end
         if output_params.PC_east
          if model.grid.Cxu < model.grid.nx
          h_out_line_e = model.fields.gh.h[model.grid.Cxu + 1,model.grid.Cyl-y_s:model.grid.Cyu+y_n]
          end
         end  
         if output_params.PC_south
          if model.grid.Cyl > 1
          h_out_line_s = model.fields.gh.h[model.grid.Cxl-x_w:model.grid.Cxu+x_e,model.grid.Cyl-1]
          end
         end
         if output_params.PC_north
          if model.grid.Cyu < model.grid.ny
          h_out_line_n = model.fields.gh.h[model.grid.Cxl-x_w:model.grid.Cxu+x_e,model.grid.Cyu+1]
          end
         end
            write_vel(simulation,h_out_line_w,h_out_line_e,h_out_line_n,h_out_line_s,x_w,x_e,y_s,y_n)
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
function write_vel(simulation::Simulation,h_out_line_w,h_out_line_e,h_out_line_n,h_out_line_s,x_w,x_e,y_s,y_n)
    @unpack model, output_params = simulation  
    
    clock_time_C=Int(round(((simulation.clock.time - simulation.timestepping_params.dt)*(3600*24*365))/output_params.dt_coup))
    clock_time_Cn=Int(round((simulation.clock.time*(3600*24*365))/output_params.dt_coup))
    
    if output_params.Output_vel    
    
    uVel_file_string = string(simulation.output_params.prefix, "_U", lpad(clock_time_C, 10,"0") ,".bin")
    vVel_file_string = string(simulation.output_params.prefix, "_V", lpad(clock_time_C, 10,"0") ,".bin")
    
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
    
    end
    
    if output_params.Output_float 
        
    float_file_string = string(simulation.output_params.prefix, "_F", lpad(clock_time_C, 10,"0") ,".bin")
    
    f_out = zeros(model.grid.Cxu - model.grid.Cxl + 1 + 2 + x_e +x_w,model.grid.Cyu - model.grid.Cyl +1 +2 + y_s +y_n)
    
    f_out[2:end-1,2:end-1]=model.fields.gh.grounded_fraction[model.grid.Cxl-x_w:model.grid.Cxu+x_e,model.grid.Cyl-y_s:model.grid.Cyu+y_n]
    
    f_out .= hton.(f_out)
    
    ffileID =  open(float_file_string,"w")
      write(ffileID, f_out[:,:])
    close(ffileID)
        
    end
    
    if output_params.Output_dhdt 
    
        dhdt_file_string = string(simulation.output_params.prefix, "_dhdt", lpad(clock_time_C, 10,"0") ,".bin")   
        
        dhdt_out = zeros(model.grid.Cxu - model.grid.Cxl + 1 + 2 + x_e +x_w,model.grid.Cyu - model.grid.Cyl +1 +2 + y_s +y_n)
        
       #doesn't include boundary rows if present
        
        dhdt_out[2+x_w:end-1-x_e,2+y_s:end-1-y_n]=model.fields.gh.dhdt[model.grid.Cxl:model.grid.Cxu,model.grid.Cyl:model.grid.Cyu]
        
        dhdt_out .= hton.(dhdt_out)
        
        ffileID =  open(dhdt_file_string,"w")
        write(ffileID, dhdt_out[:,:])
        close(ffileID)
        
    end
   
    if output_params.Output_Hb
    
    if model.grid.Cxl > 1 || model.grid.Cxu <  model.grid.nx || model.grid.Cyl > 1 || model.grid.Cyu < model.grid.ny 
            
            
     h_out_b = zeros(model.grid.Cxu - model.grid.Cxl + 1 + 2 + x_e +x_w ,model.grid.Cyu - model.grid.Cyl +1 + 2 + y_s + y_n)
     if output_params.PC_west       
      if model.grid.Cxl > 1
                h_out_b[2,2:end-1] .= h_out_line_w[:]
      end
     end
     if output_params.PC_east
      if model.grid.Cxu < model.grid.nx
                h_out_b[end-1,2:end-1] .= h_out_line_e[:]
      end
     end
     if output_params.PC_south
      if model.grid.Cyl > 1
                h_out_b[2:end-1,2] .= h_out_line_s[:]
      end
     end
     if output_params.PC_north
      if model.grid.Cyu < model.grid.ny
                h_out_b[2:end-1,end-1] .= h_out_line_n[:]
      end
     end
     h_out_b .= hton.(h_out_b)
     
        
     hb_file_string = string(simulation.output_params.prefix, "_Hb", lpad(clock_time_Cn, 10,"0") ,".bin")
    
     hbfileID =  open(hb_file_string,"w")
     write(hbfileID, h_out_b[:,:])
     close(hbfileID) 
     
     end
        
    end
 end 
