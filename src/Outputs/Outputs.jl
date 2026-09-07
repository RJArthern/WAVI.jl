module Outputs

export OutputParams, with_cleared_stencil_scratch, is_output_step

import WAVI.Deferred: clear!, collect!
using WAVI.Deferred
using Parameters
using WAVI: AbstractSpec
using WAVI.Time: Clock

#structure that contains outputting info
struct OutputParams{T<:Real, R<:Real, O<:Collector}
    outputs::O                  # Data Collection of outputs that can be lazily output
    output_freq::T              # output time 
    n_iter_out::R               # number of steps per output
    output_format::String       # specify output format [mat/jld]
    prefix::String              # file prefix
    output_path::String         # folder in which to save
    dump_vel::Bool              # toggle on dumping the velocity after the final timestep
    zip_format::String          # specify whether or not to zip the output, and the format
    output_start::Bool          # flag to specify whether to output the initial state or not
end

"""
    OutputParams(; 
        outputs = (),
        output_freq = Inf, 
        output_format = "jld2",
        prefix = "outfile", 
        output_path = "./",
        dump_vel = false,
        zip_format = "none",
        output_start = false)

Construct a WAVI.jl output parameters object.

Keyword arguments
=================
- `output_freq`: quantity specifying hwo frequently to produce output 
- `output_format`: specify output format (currently only .mat and .jld2 outputs are supported, selected with 'mat' or 'jld2' options)
- `prefix`: prefix to be prepended onto file names
- `output_path`: path on which to output filese
- `dump_vel`: flag to toggle whether or not dumping the velocity after the final timestep (used in mitgcm coupling)
- `zip_format`: specify whether or not to zip the output, and the format (currently only '.nc' output supported, use flag 'nc')
- `output_start`: flag to specify whether to output at the zeroth time step    
"""
function OutputParams(outputs::NamedTuple; 
    output_freq = Inf, 
    output_format = "jld2",
    prefix = "outfile", 
    output_path = "./",
    dump_vel = false,
    zip_format = "none",
    output_start = false)

    #default the n_iter_out to -1 (this is updated in simulation once we know timestep from timestepping_params)
    n_iter_out = -1

    #check output_freq
    ((output_freq == Inf) || (output_freq > 0)) || throw(ArgumentError("output frequency must be positive or Inf"))

    #create output path if it doesn't exist (using mkpath instead of mkdir to avoid MPI races)
    mpi_rank = get(ENV, "OMPI_COMM_WORLD_RANK", get(ENV, "PMI_RANK", "0"))
    is_root_rank = mpi_rank == "0"
    if !isdir(output_path) && is_root_rank
        @warn string("Did not find output path ", output_path, ", creating it")
    end
    mkpath(output_path)

    #append a "/" to folder if it doesn't have one
    endswith(output_path, "/") || (output_path = string(output_path, "/"))

    #throw an error if we don't get an output format we know
    ((output_format == "jld2") || (output_format == "mat")) || throw(ArgumentError("Output format must be `jld2` or `mat`"))

    #revert the zip format to none if we don't recognise format
    if ~(zip_format in ["none", "nc"])
        println("detected a zip format other than none or nc...
        WAVI currently only supports zipping to nc.
        Reverting to no zipping")
        zip_format = "none"
    end

    collector = Collector()
    for (nom, path) in pairs(outputs)
        register_field!(collector, nom, path)
    end
    return OutputParams(collector, output_freq, n_iter_out, output_format, prefix, output_path, dump_vel, zip_format, output_start)
end

function OutputParams(; outputs = NamedTuple(), kwargs...)
    return OutputParams(outputs; kwargs...)
end
clear!(op::OutputParams) = clear!(op.outputs)
collect!(op::OutputParams, args...) = collect!(op.outputs, args...)

"""
    is_output_step(output_params, clock)

True if this step should write field output.

That is when `n_iter_out` (steps between writes) is positive and `clock.n_iter` is
a multiple of it.
"""
function is_output_step(output_params::OutputParams, clock::Clock)
    n = output_params.n_iter_out
    return n > 0 && mod(clock.n_iter, n) == 0
end

include("checkpoints.jl")
include("output_writing.jl")
include("zipping_output.jl")

end