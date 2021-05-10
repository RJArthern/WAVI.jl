
#structure that contains outputting info
struct OutputParams{T <: Real, R <: Real, O}
    outputs::O         #tuple of entries defining names and quantities to be outputted
    output_freq::T     #output time 
    n_iter_out::R      #number of steps per output
    format::String     #specify output format [mat/jld]
    zipped::Bool       #dump files at each timestep or zip the files
    prefix::String     #file prefix
end

#output constructor
function OutputParams(; 
    outputs = (),
    output_freq = Inf, 
    format = "jld2",
    zipped = false,
    prefix = "outfile")

    #default the n_iter_out to Inf (this will be updated in the simulation)
    n_iter_out = Inf

    #check the format
    (format == "jld2" || format == "mat") || ArgumentError("Output format must be jld2 or mat")

    return OutputParams(outputs, output_freq, n_iter_out, format, zipped, prefix)
end
