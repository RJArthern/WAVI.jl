#file containing outputting functions
"""
    write_output(simulation::AbstractSimulation)

Output the data from the simulation at the current timestep
"""
function write_output(simulation::AbstractSimulation)
    @unpack output_params, model, clock = simulation
    output_dict = fetch_output(output_params.outputs)

    #put the grid co-ordinates and time into output
    if ~haskey(output_dict, :t); output_dict["t"] = clock.time; end
    if ~haskey(output_dict, :x); output_dict["x"] = model.grid.xxh; end
    if ~haskey(output_dict, :y); output_dict["y"] = model.grid.yyh; end

    fname = string(output_params.prefix , lpad(simulation.clock.n_iter, 10,"0"));
    if output_params.format == "jld2"
        fname = string(fname, ".jld2")
        save(fname, output_dict)
    elseif output_params.format == "mat"
        fname = string(fname, ".mat")
        matwrite(fname, output_dict)
    end
end

"""
    fetch_output(outputs)

Return a dictionary with dictionary entries corresponding to outputs
"""
function fetch_output(outputs)
    output_dict = Dict()
    for (k,v) in zip(keys(outputs), outputs)
        output_dict[string(k)] = v
    end
    return output_dict
end

