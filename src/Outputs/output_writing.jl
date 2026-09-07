
export write_outputs, write_output, write_vel

using JLD2
using MAT

using WAVI
using WAVI: AbstractModel, AbstractSimulation
using WAVI.Deferred
using WAVI.Parameters
using WAVI.Time

function write_outputs(model::M,
                       timestepping_params::TimesteppingParams, 
                       output_params::OutputParams, 
                       clock::Clock) where {M<:AbstractModel{<:Any, <:Any, <:Any}}

    @unpack spec = model
                       
    if should_write_checkpoint(timestepping_params, clock)
        write_checkpoint!(spec, model, timestepping_params, output_params, clock)
    end

    #check if we have hit an output timestep
    if is_output_step(output_params, clock)
        write_output(model, output_params, clock)
    end

    #check the dump velocity flag at the final timestep
    if (clock.n_iter == timestepping_params.n_iter_total) && output_params.dump_vel
        write_vel(output_params, model)
    end
end
write_outputs(sim::AbstractSimulation) = write_outputs(sim.model, sim.timestepping_params, sim.output_params, sim.clock)

#file containing outputting functions
"""
    write_output(model, output_params, clock)

Output the data from the simulation at the current timestep
"""
function write_output(model::M, output_params::OutputParams, clock::Clock) where {M<:AbstractModel{<:Any, <:Any, <:Any}}
    output_dict = collect!(output_params, model)
    name = lpad(clock.n_iter, 10,"0")

    if isnothing(output_dict)
        @warn "No outputs processed for $(name)"
    end

    #put the grid co-ordinates and time into output.
    #Round time in output to some decimal places to make it prettier (machine precision can make this look nasty!)
    if ~haskey(output_dict, :t); output_dict["t"] = round(clock.time, digits = 3); end
    if ~haskey(output_dict, :x); output_dict["x"] = model.grid.xxh; end
    if ~haskey(output_dict, :y); output_dict["y"] = model.grid.yyh; end
    _host_copy_output!(output_dict)

    fname = string(output_params.output_path, output_params.prefix, name)
    if output_params.output_format == "jld2"
        fname = string(fname, ".jld2")
        save(fname, output_dict)
    elseif output_params.output_format == "mat"
        fname = string(fname, ".mat")
        matwrite(fname, output_dict)
    end
    
    @info "Output at timestep number $(clock.n_iter) - $(fname)"
    clear!(output_params)
end

# NetCDF and MAT writers index arrays on the host. Copy device fields once here
# so JLD2 dumps are also CPU-readable found by error on me testing outputs.
function _host_copy_output!(d::AbstractDict)
    for (k, v) in d
        if v isa AbstractArray && !(v isa Array || v isa BitArray)
            d[k] = Array(v)
        end
    end
    return d
end

"""
    function write_vel(simulation)

Write the velocity at the the final timestep of the simulation (used in the coupled wavi-mitgcm model to communicate with streamice)

    TODO: This can be made more generic for binary dumps of any field
"""
function write_vel(output_params::OutputParams, model::AbstractModel)
    uVel_file_string = string(output_params.prefix,  "_U.bin")
    vVel_file_string = string(output_params.prefix,  "_V.bin")
    
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


