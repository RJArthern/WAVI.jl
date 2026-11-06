export get_spatiotemporal_var_atts, get_spatial_dimensions, get_times, get_output_as_dict, 
    make_ncfile, make_ncfile_from_filenames, zip_output

using MAT
using NCDatasets

using WAVI: AbstractModel, AbstractSimulation
using WAVI.Outputs: OutputParams

"""
    get_format_filenames(format, folder)

Returns an array of filenames in folder (string) with suffix format (string)
"""
get_format_filenames(format::String, folder::String, prefix::String) =[string(folder,f) for f in  readdir(folder) if (endswith(f, format) && startswith(f,prefix))]


"""
    get_spatiotemporal_var_atts()

Return the variable attributes for the spatiotemporal variable (x,y,y)
"""
function get_spatiotemporal_var_atts()
    x_atts = Dict("longname" => "x co-ordinates of ice grid points (h grid)",  "units" => "m")
    y_atts = Dict("longname" => "y co-ordinates of ice grid points (h grid)",  "units" => "m")
    time_atts = Dict("longname" => "Time", "units" => "years");
    return x_atts, y_atts, time_atts
end

"""
    get_spatial_dimensions()

Return one-dimensional arrays of the spatial variables
"""
function get_spatial_dimensions(fname)
    format = return_extension(fname)
    if format == "mat"
        vars = matread(fname)
    elseif format == "jld2"
        vars = load(fname)
    end
    X = vars["x"]
    X =  X[:,1]
    Y = vars["y"]
    Y = Y[1,:]
    return X, Y
end

"""
    get_time_output(filenames)

Return the times associated with filenames
"""

function get_times(filenames)
    t = zeros(length(filenames))
    for i = 1:length(filenames)
        format = return_extension(filenames[i])
        vars = get_output_as_dict(filenames[i],format)
        t[i] = vars["t"]       
    end
    return t
end
"""
    get_output_as_dict(filename,format)

Read the data in filename according to different format specification
"""
function get_output_as_dict(filename,format)
    if format == "mat"
        vars= matread(filename)
    elseif format == "jld2"
        vars = load(filename)
    end
return vars
end

"""
    return_extension(file)

Return the extension of the input file 
"""
return_extension(file) =  file[(findlast(isequal('.'),file)+1):end];

"""
make_ncfile(folder,format)

Wrapper script to zip the output files in "folder" with type "format" to an nc file with name nc_name_full (including path)
"""
function make_ncfile(format, folder, nc_name, prefix)
    #check that the input format is 
    filenames = get_format_filenames(format, folder, prefix)
    if ~isempty(filenames)
        make_ncfile_from_filenames(filenames, format, nc_name)
    else
        println("attempted to zip the outputs to nc format, but did not find any files...")
    end
    return nothing
end

"""
    make_ncfile_from_filenames(filenames, format)

Output an nc file from filenames, which have format "format" to a file with name nc_name.
nc_name must contain the path as well!
"""
function make_ncfile_from_filenames(filenames, format, nc_name_full)
    #get the spatial and temporal variables as arrays of size N x 1
    x, y = get_spatial_dimensions(filenames[1])
    t    = get_times(filenames)

    #setup attributes for spatiotemporal variables 
    x_atts, y_atts, time_atts = get_spatiotemporal_var_atts()

    # Load all variables from first snapshot file to inspect their
    # sizes and types without reading files more than once
    first_file_vars = get_output_as_dict(filenames[1], format)

    # Keep only data variables, excluding co-ordinate keys (x, y, t)
    data_keys = filter(k -> !(k in ["x", "y", "t"]), collect(keys(first_file_vars)))

    # Get the size of sigma (vertical) dimension, if any 3D fields are present
    # Sigma not stored explicitly in snapshot files, so infer from shape of first 3D variable
    nσ = nothing
    for k in data_keys
        if ndims(first_file_vars[k]) == 3
            nσ = size(first_file_vars[k], 3)
            break
        end
    end

    # Remove existing file if it exists
    isfile(nc_name_full) && rm(nc_name_full)
    
    NCDataset(nc_name_full, "c") do ds
        defDim(ds, "x", length(x))
        defDim(ds, "y", length(y))
        defDim(ds, "TIME", length(t))
        !isnothing(nσ) && defDim(ds, "σ", nσ)

        # Define and write coordinate variables
        defVar(ds, "x", x, ("x",), attrib = x_atts)
        defVar(ds, "y", y, ("y",), attrib = y_atts)
        defVar(ds, "TIME", t, ("TIME",), attrib = time_atts)

        for key in data_keys
            # Use the already-loaded first-file data to check the variable's type and size.
            first_file_data = first_file_vars[key]
            sz = size(first_file_data)

            # NetCDF writer does not support Bool directly; store logical fields as UInt8 (0/1).
            is_bool = eltype(first_file_data) <: Bool
            var_type = is_bool ? UInt8 : eltype(first_file_data)

            if sz == (length(x), length(y))
                # 2D field: dimensions are (x, y, TIME)
                var_nc = defVar(ds, key, var_type, ("x", "y", "TIME"))
                # Populate the variable by iterating through filenames
                for i in 1:length(filenames)
                    data_i = get_output_as_dict(filenames[i], format)[key]
                    var_nc[:, :, i] = is_bool ? UInt8.(data_i) : data_i
                end
            elseif !isnothing(nσ) && sz == (length(x), length(y), nσ)
                # 3D field: dimensions are (x, y, sigma, TIME)
                var_nc = defVar(ds, key, var_type, ("x", "y", "σ", "TIME"))
                for i in 1:length(filenames)
                    data_i = get_output_as_dict(filenames[i], format)[key]
                    var_nc[:, :, :, i] = is_bool ? UInt8.(data_i) : data_i
                end
            else
                @warn string("found an output variable (", key, ") who's spatial dimensions (", sz, ") do not match the co-ordinates. Skipping this variable from the nc output...")
            end
        end
    end
    return nothing
end


"""
    zip_output(simulation)

Zip all of the output files from simulation.
"""
function zip_output(model::M, output_params::OutputParams) where {M<:AbstractModel{<:Any, <:Any, <:Any}}
    if output_params.zip_format == "nc"
        nc_name_full = string(output_params.output_path, output_params.prefix, ".nc")
        @info "Creating NetCDF output $(nc_name_full)"
        make_ncfile(output_params.output_format, output_params.output_path, nc_name_full, output_params.prefix)
    end
    return nothing
end
zip_output(sim::AbstractSimulation) = zip_output(sim.model, sim.output_params)
