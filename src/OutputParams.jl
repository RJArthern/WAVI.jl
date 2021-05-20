
#structure that contains outputting info
struct OutputParams{T <: Real, R <: Real, O}
    outputs::O         #tuple of entries defining names and quantities to be outputted
    output_freq::T     #output time 
    n_iter_out::R      #number of steps per output
    format::String     #specify output format [mat/jld]
    prefix::String     #file prefix
    output_path::String     #folder in which to save
    dump_vel::Bool     #toggle on dumping the velocity after the final timestep
end

#output constructor
function OutputParams(; 
    outputs = (),
    output_freq = Inf, 
    format = "jld2",
    prefix = "outfile", 
    output_path = "./",
    dump_vel = false)

    #default the n_iter_out to Inf (this will be updated in the simulation)
    n_iter_out = Inf

    #if you don't find folder, set it to the working directory
    if ~isdir(output_path)
        @warn string("Did not find output path ", output_path, ". Any outputs will go to the working directory", pwd())
        output_path = "./"
    end

    #append a "/" to folder if it doesn't have one
    endswith(output_path, "/") || (output_path = string(output_path, "/"))

    #check the format
    (format == "jld2" || format == "mat") || ArgumentError("Output format must be jld2 or mat")

    return OutputParams(outputs, output_freq, n_iter_out, format, prefix, output_path, dump_vel)
end
