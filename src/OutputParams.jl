
#structure that contains outputting info
struct OutputParams{T <: Real, R <: Real}
    out_dict::Dict       #dictionary of entries defining names and quantities to be outputted
    out_dict_values::Dict #dictionary that stores the evaluated form of the output dictionary
    out_freq::T         #output time 
    n_iter_out::R       #number of steps per output
    format::String      #specify output format [mat/jld]
    zipped::Bool        #dump files at each timestep or zip the files
    prefix::String      #file prefix
end

#output constructor
function OutputParams(; 
    out_dict = Dict(),
    out_freq = Inf, 
    format = "jld2",
    zipped = false,
    prefix = "outfile")

    #default the n_iter_out to Inf (this will be updated in the simulation)
    n_iter_out = Inf

    #default the output dict values to empty
    out_dict_values = Dict()

    #check the format
    (format == "jld2" || format == "mat") || ArgumentError("Output format must be jld2 or mat")

    return OutputParams(out_dict, out_dict_values,out_freq, n_iter_out, format, zipped, prefix)
end
