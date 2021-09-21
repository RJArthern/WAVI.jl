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

    #setup Dimensions
    t_ncdim = NcDim("TIME", t, atts = time_atts)
    x_ncdim = NcDim("x", x, atts = x_atts)
    y_ncdim = NcDim("y", y, atts = y_atts)

    #create the nc file variables
    output_dict = Dict() #dictionary for the output values
    output_vars = Dict() #dictionary for the output variables
    
    #get the keys
    if format == "mat"
        filekeys = keys(matread(filenames[1]))
    elseif format == "jld2"
        filekeys = keys(load(filenames[1]))
        #println(typeof(filekeys))
    end

    for key in filekeys
        if ~(key in ["x", "y", "t"]) #if this isn't a spatial dimensions
            var_out = zeros(length(x), length(y), length(t))
    
            #check the size of this variable in the first file
            sz = size(get_output_as_dict(filenames[1],format)[key])
            if sz == (length(x), length(y))
                #if the size is OK, create a variable and add to the array
                ncvar_key = NcVar(key, Array([x_ncdim, y_ncdim, t_ncdim]))
                output_vars[key] = ncvar_key
                
                #get the data in an array
                for i = 1:length(filenames)
                    var_out[:,:,i] = get_output_as_dict(filenames[i],format)[key]
                end

                #add this to dictionary
                output_dict[key] = var_out
            else
                @warn string("found an output variable (", key, ") who's spatial dimensions do not match the co-ordinates. Skipping this variable from the nc output...")
            end
        end
    end

    #make the nc file
    isfile(nc_name_full) && rm(nc_name_full)
    NetCDF.create(nc_name_full, [val for val in values(output_vars)]) do nc
        # Writing data to the file is done using putvar
        for key in keys(output_dict) #for every variable of the correct size
            NetCDF.putvar(nc,key,output_dict[key])
        end
    end
    return nothing
end


"""
    zip_output(simulation)

Zip all of the output files from simulation.
"""
function zip_output(simulation)
    @unpack output_params = simulation
    if output_params.zip_format == "nc"
        nc_name_full = string(output_params.output_path, output_params.prefix, ".nc")
        make_ncfile(output_params.output_format, output_params.output_path, nc_name_full, output_params.prefix)
    end
    return nothing
end