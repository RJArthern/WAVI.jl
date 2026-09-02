"""
    get_format_filenames(format, folder)

Returns an array of filenames in folder (string) with suffix format (string)
"""
get_format_filenames(format::String, folder::String, prefix::String) =[string(folder,f) for f in  readdir(folder) if (endswith(f, format) && startswith(f,prefix))]


"""
    detect_grid_type(dims, grid_dims)

Detects the grid type (h, u, v) based on variable dimensions.
Returns the grid type as a symbol or :unknown if no match.
"""
function detect_grid_type(dims, grid_dims)
    if dims == grid_dims[:h]
        return :h
    elseif dims == grid_dims[:u]
        return :u
    elseif dims == grid_dims[:v]
        return :v
    else
        return :unknown
    end
end

"""
    get_grid_dimensions(fname)

Extract grid dimensions for h, u, and v grids from the first output file.
Returns a dictionary with grid dimensions for each grid type.
"""
function get_grid_dimensions(fname)
    format = return_extension(fname)
    if format == "mat"
        vars = matread(fname)
    elseif format == "jld2"
        vars = load(fname)
    end
    
    # Extract dimensions from coordinate arrays
    nx, ny = size(vars["x"])
    
    grid_dims = Dict(
        :h => (nx, ny),
        :u => (nx + 1, ny),
        :v => (nx, ny + 1)
    )
    
    return grid_dims
end

"""
    get_spatiotemporal_var_atts(grid_type::Symbol)

Return the variable attributes for the spatiotemporal variable (x,y,t) for a specific grid type.
"""
function get_spatiotemporal_var_atts(grid_type::Symbol)
    grid_name = string(grid_type)
    x_atts = Dict("long_name" => "x coordinates of grid points ($grid_name grid)", "units" => "m")
    y_atts = Dict("long_name" => "y coordinates of grid points ($grid_name grid)", "units" => "m")
    return x_atts, y_atts
end

"""
    get_spatial_dimensions_for_grid(fname, grid_type::Symbol)

Return one-dimensional arrays of the spatial variables for a specific grid type.
"""
function get_spatial_dimensions_for_grid(fname, grid_type::Symbol)
    format = return_extension(fname)
    if format == "mat"
        vars = matread(fname)
    elseif format == "jld2"
        vars = load(fname)
    end
    
    if grid_type == :h
        X = vars["x"][:,1]
        Y = vars["y"][1,:]
    elseif grid_type == :u
        # Check if u-grid coordinates exist, otherwise create them from h-grid
        if haskey(vars, "xu") && haskey(vars, "yu")
            X = vars["xu"][:,1]
            Y = vars["yu"][1,:]
        else
            # Create u-grid coordinates from h-grid coordinates
            xh = vars["x"][:,1]
            yh = vars["y"][1,:]
            dx = length(xh) > 1 ? xh[2] - xh[1] : 8000.0
            X = [xh[1] - dx/2; xh .+ dx/2]  # u-grid has nx+1 points
            Y = yh  # y coordinates same as h-grid
        end
    elseif grid_type == :v
        # Check if v-grid coordinates exist, otherwise create them from h-grid
        if haskey(vars, "xv") && haskey(vars, "yv")
            X = vars["xv"][:,1]
            Y = vars["yv"][1,:]
        else
            # Create v-grid coordinates from h-grid coordinates
            xh = vars["x"][:,1]
            yh = vars["y"][1,:]
            dy = length(yh) > 1 ? yh[2] - yh[1] : 8000.0 # default spacing 
            X = xh  # x coordinates same as h-grid
            Y = [yh[1] - dy/2; yh .+ dy/2]  # v-grid has ny+1 points
        end
    else
        throw(ArgumentError("Unknown grid type: $grid_type"))
    end
    
    return X, Y
end

"""
    get_times(filenames)

Extract time values from the output files.
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
nc_name must contain the path as well! Supports multiple grids (h, u, v) with grid type detection.
"""
function make_ncfile_from_filenames(filenames, format, nc_name_full)
    # Get time
    t = get_times(filenames)
    
    # Get the grid dimensions for type detection
    grid_dims = get_grid_dimensions(filenames[1])

    # Get the keys of variables to be written (excluding spatial and time dimensions)
    filekeys = collect(keys(get_output_as_dict(filenames[1], format)))
    data_keys = filter(k -> !(k in ["x", "y", "xu", "yu", "xv", "yv", "t"]), filekeys)
    
    # Determine which grids are actually used by analyzing variable dimensions
    used_grids = Set{Symbol}()
    first_file_vars = get_output_as_dict(filenames[1], format)
    
    for key in data_keys
        dims = size(first_file_vars[key])
        grid_type = detect_grid_type(dims, grid_dims)
        if grid_type != :unknown
            push!(used_grids, grid_type)
        end
    end
    
    # Only extract coordinates for grids that are actually used
    grid_coords = Dict{Symbol, Tuple{Vector, Vector}}()
    for grid_type in used_grids
        try
            X, Y = get_spatial_dimensions_for_grid(filenames[1], grid_type)
            grid_coords[grid_type] = (X, Y)
        catch e
            @warn "Could not extract coordinates for $grid_type grid: $e"
        end
    end

    # Remove existing file if it exists
    isfile(nc_name_full) && rm(nc_name_full)
    
    NCDataset(nc_name_full, "c") do ds
        # Define time dimension
        defDim(ds, "TIME", length(t))
        time_atts = Dict("long_name" => "Time", "units" => "years")
        defVar(ds, "TIME", t, ("TIME",), attrib = time_atts)
        
        # Define dimensions and coordinate variables for each available grid
        for (grid_type, (X, Y)) in grid_coords
            x_dim_name = "x_$(grid_type)"
            y_dim_name = "y_$(grid_type)"
            x_var_name = "x_$(grid_type)"
            y_var_name = "y_$(grid_type)"
            
            # For backward compatibility, keep the original dimension and variable names for h-grid
            if grid_type == :h
                x_dim_name = "x"
                y_dim_name = "y"
                x_var_name = "x"
                y_var_name = "y"
            end
            
            # Define dimensions
            defDim(ds, x_dim_name, length(X))
            defDim(ds, y_dim_name, length(Y))
            
            # Get attributes for this grid type
            x_atts, y_atts = get_spatiotemporal_var_atts(grid_type)
            
            # Define coordinate variables
            defVar(ds, x_var_name, X, (x_dim_name,), attrib = x_atts)
            defVar(ds, y_var_name, Y, (y_dim_name,), attrib = y_atts)
        end

        # Process each data variable
        for key in data_keys
            # Get the data from the first file to determine its type and size
            first_file_data = get_output_as_dict(filenames[1], format)[key]
            sz = size(first_file_data)
            
            # Detect which grid this variable belongs to
            grid_type = detect_grid_type(sz, grid_dims)
            
            if grid_type != :unknown && haskey(grid_coords, grid_type)
                # Convert Bool type to Int8 for NetCDF compatibility
                data_type = eltype(first_file_data)
                if data_type == Bool
                    data_type = Int8
                end
                
                # Use appropriate dimension names for this grid
                x_dim_name = "x_$(grid_type)"
                y_dim_name = "y_$(grid_type)"
                
                # For backward compatibility, use original dimension names for h-grid
                if grid_type == :h
                    x_dim_name = "x"
                    y_dim_name = "y"
                end
                
                # Define the variable in the NetCDF file
                var_nc = defVar(ds, key, data_type, (x_dim_name, y_dim_name, "TIME"))
                
                # Add grid type as an attribute
                var_nc.attrib["grid_type"] = string(grid_type)

                # Populate the variable by iterating through filenames
                for i = 1:length(filenames)
                    data = get_output_as_dict(filenames[i], format)[key]
                    # Convert Bool arrays to Int8 for NetCDF compatibility
                    if eltype(data) == Bool
                        var_nc[:,:,i] = Int8.(data)
                    else
                        var_nc[:,:,i] = data
                    end
                end
                
                println("Added variable '$key' to NetCDF file on $(grid_type) grid with dimensions $sz")
            else
                @warn string("found an output variable (", key, ") with spatial dimensions (", sz, ") that do not match any known grid type (h: $(grid_dims[:h]), u: $(grid_dims[:u]), v: $(grid_dims[:v])). Skipping this variable from the nc output...")
            end
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