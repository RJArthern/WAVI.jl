#file containing outputting functions

#structure that contains outputting info
struct Output{T <: Real, N}
    out_dict::Dict       #dictionary of entries defining names and quantities to be outputted
    out_freq::T         #output time 
    n_iter_out::N       #number of steps per output
    format::String      #specify output format [mat/jld]
    zipped::Bool        #dump files at each timestep or zip the files
    prefix::String      #file prefix
end

#output constructor
function Output(; 
    out_dict = Dict(),
    out_freq = Inf, 
    format = "jld2",
    zipped = false,
    prefix = "outfile")

    #default the n_iter_out to Inf (this will be updated in the simulation)
    n_iter_out = Inf

    return Output(out_dict, out_freq, n_iter_out, format, zipped, prefix)
end

#function write_output(wavi::AbstractModel,output::Output