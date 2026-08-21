using JLD2
using MAT
using MPI

using WAVI.Parameters

import WAVI: AbstractModel
import WAVI.Grids: Grid
import WAVI.MeltRates: UniformMeltRate
import WAVI.Models: BasicSpec, Model, get_bed_elevation
import WAVI.Outputs: write_outputs, zip_output, OutputParams, should_write_checkpoint, write_checkpoint!
import WAVI.Parameters: TimesteppingParams
import WAVI.Time: Clock


# TODO: these redefinitions are loathsome, but can't get a clearer way of limiting output via a simple method declaration with @root
function write_output(model::M, output_params::OutputParams, clock::Clock) where {M<:AbstractModel{<:Any, <:Any, <:MPISpec}}
    output_dict = collect!(output_params, model)
    name = lpad(clock.n_iter, 10,"0")

    @root begin
        if isnothing(output_dict)
            @warn "No outputs processed for $(name)"
        end

        #put the grid co-ordinates and time into output.
        #Round time in output to some decimal places to make it prettier (machine precision can make this look nasty!)
        if ~haskey(output_dict, :t); output_dict["t"] = round(clock.time, digits = 3); end
        if ~haskey(output_dict, :x); output_dict["x"] = model.global_grid.xxh; end
        if ~haskey(output_dict, :y); output_dict["y"] = model.global_grid.yyh; end

        fname = string(output_params.output_path, output_params.prefix, name)
        if output_params.output_format == "jld2"
            fname = string(fname, ".jld2")
            save(fname, output_dict)
        elseif output_params.output_format == "mat"
            fname = string(fname, ".mat")
            matwrite(fname, output_dict)
        end
        
        @info "Output at timestep number $(clock.n_iter) - $(fname)"
    end
    clear!(output_params)
end

function write_outputs(model::M,
                       timestepping_params::TimesteppingParams, 
                       output_params::OutputParams, 
                       clock::Clock) where {M<:AbstractModel{<:Any, <:Any, <:MPISpec}}

    if should_write_checkpoint(timestepping_params, clock)
        write_checkpoint!(model, timestepping_params, output_params, clock)
    end

    #check if we have hit an output timestep
    if mod(clock.n_iter, output_params.n_iter_out) == 0
        write_output(model, output_params, clock)
    end

        #check the dump velocity flag at the final timestep
#         if (clock.n_iter == timestepping_params.n_iter_total) && output_params.dump_vel
#             write_vel(output_params, model)
#         end
end

function zip_output(model::M, output_params::OutputParams) where {M<:AbstractModel{<:Any, <:Any, <:MPISpec}}
    @info "[$(model.spec.rank + 1)/$(model.spec.global_size)] Called zip output"
    @root begin
        if output_params.zip_format == "nc"
            nc_name_full = string(output_params.output_path, output_params.prefix, ".nc")
            @info "Creating NetCDF output on MPI root $(nc_name_full)"
            make_ncfile(output_params.output_format, output_params.output_path, nc_name_full, output_params.prefix)
        end
    end
    return nothing
end
