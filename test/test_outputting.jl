using Test, WAVI, MAT, JLD2, NCDatasets, Dates

function output_test(; dt  = 0.5, 
                    end_time = 100., 
                    output_freq = 5., 
                    output_format = "jld2", 
                    zip_format = "none", 
                    prefix = "outfile", 
                    dump_vel = false, 
                    output_path = ".",
                    pchkpt_freq = Inf,
                    chkpt_freq = Inf,
                    niter0 = 0,
                    output_start = false)

    grid = Grid() #default grid with nx = 80, ny = 10
    bed = WAVI.mismip_plus_bed
    solver_params = SolverParams(maxiter_picard = 1)
    model = Model(grid = grid, bed_elevation = bed, solver_params = solver_params)
    timestepping_params = TimesteppingParams(niter0 = niter0, dt = dt, end_time = end_time, chkpt_freq = chkpt_freq, pchkpt_freq = pchkpt_freq)

    outputs = (h = model.fields.gh.h,
                u = model.fields.gh.u,
                v = model.fields.gh.v,
                b = model.fields.gh.b) #output velocities and thickness

    output_freq = output_freq
    output_params = OutputParams(outputs = outputs, 
                        output_freq = output_freq,
                        output_format = output_format,
                        output_path = output_path,
                        zip_format = zip_format,
                        prefix = prefix,
                        dump_vel = dump_vel,
                        output_start = output_start)

    simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params,
                        output_params = output_params)
            
    #perform the simulation
    run_simulation!(simulation)
    return simulation
end

#########################################################
############### Test flags ##############################
test_basic_outputting = true
test_zipping_output   = true
test_output_after_restart = true
test_output_errors = true
test_multigrid_output = true
##########################################################


@testset "Outputting" begin
    if test_basic_outputting
        @testset "Output files" begin 
            @info "Testing outputting..."
            
            for output_format in ["mat", "jld2"] #check that both mat and jld2 work
                for folder = ["outputs/", "outputs", "./outputs"] #check that both with and without / works
                    isdir(folder) && rm(folder, force = true, recursive = true) #make a new folder
                    mkdir(folder)

                    #run the simulation
                    sim = output_test(output_path = folder, 
                                end_time = 100., 
                                output_freq = 5., 
                                prefix = "testoutfile", 
                                output_format = output_format)

                    foldersim = sim.output_params.output_path 
                    @test foldersim[end] == '/' #test that we do have the / at end of path
                    files = [string(foldersim,f) for f in  readdir(foldersim) if endswith(f, output_format)]
                    println(files)
                    @test length(files) == 20    #check there are the correct number of output files

                    #check they have the correct suffix
                    prefices = [split(str,".")[1] for str in files]
                    suffices = [split(str,".")[end] for str in files]
                    @test all(suffices .== output_format)

                    #check that variables included and correct size
                    fname = string(foldersim,readdir(folder)[1])
                    if output_format == "mat"
                    dict = matread(fname)
                    elseif output_format == "jld2"
                    dict = load(fname)
                    end
                    @test length(dict["t"]) == 1
                    @test size(dict["x"]) == (80,10)
                    @test size(dict["y"]) == (80,10)
                    @test size(dict["b"]) == (80,10)
                    @test size(dict["u"]) == (80,10)
                    @test size(dict["h"]) == (80,10)
                    @test size(dict["v"]) == (80,10)

                    #check it has the correct prefix and number of zeros
                    fname = readdir(folder)[1]
                    @test fname[1:11]  == "testoutfile"

                    #delete the folder
                    rm(folder, force = true, recursive = true)
                end
            end #end loop over output format
        
        end
        @testset "Test output start" begin 
            for output_start in [true, false]
                folder = "outputs/"
                isdir(folder) && rm(folder, force = true, recursive = true) #make a new folder
                mkdir(folder)
                #test that if we only get an output at the start (e.g. if the output freq > end time), then outputted solution matches the expected thing
                sim = output_test(output_path = folder, 
                                end_time = 1., 
                                output_freq = 5., 
                                prefix = "testoutfile", 
                                output_format = "jld2", 
                                output_start = output_start)

                foldersim = sim.output_params.output_path 
                if output_start
                    #check only a single file outputted
                    @test length(readdir(foldersim)) == 1

                    fname = string(foldersim,readdir(folder)[1])
                    dict = load(fname)

                    #check that thickness same as IC and time = 0
                    @test dict["h"] == sim.model.initial_conditions.initial_thickness
                    @test dict["t"] == 0.0


                else
                    #check that no file output
                    @test isempty(readdir(foldersim))
                end

                #delete the folder
                rm(folder, force = true, recursive = true)
            end
        end
    end

    if test_zipping_output
        @testset "Zipping output" begin 
            @info "Testing zipping output..."
            for output_format in ["mat", "jld2"]
                for output_start in [true, false]
                    folder = "outputs/"
                    isdir(folder) && rm(folder, force = true, recursive = true) #make a new folder
                    mkdir(folder)

                    sim = output_test(zip_format = "nc", 
                                        output_path = folder, 
                                        output_format = "mat",
                                        output_start = output_start)
                    
                    @test sim isa Simulation
                    fname = string(folder, sim.output_params.prefix, ".nc")
                    @test isfile(fname) #check the zipped file exists

                    #test variables read from nc file
                    NCDataset(fname, "r") do ds
                    # Read all necessary variables from the open dataset
                        x = ds["x"][:]
                        y = ds["y"][:]
                        h = ds["h"][:,:,:]
                        u = ds["u"][:,:,:]
                        v = ds["v"][:,:,:]
                        b = ds["b"][:,:,:]
                        t = ds["TIME"][:]

                        if output_start
                            @test t == 0.:5.:100.
                        else
                            @test t == 5.:5.:100.
                        end
                            
                        @test size(x) == (80,)
                        @test size(y) == (10,)
                        @test size(b) == (length(x), length(y), length(t))
                        @test size(u) == (length(x), length(y), length(t))
                        @test size(v) == (length(x), length(y), length(t))
                        @test size(h) == (length(x), length(y), length(t))
                    end
            
                    #delete the folder
                    rm(folder, force = true, recursive = true)
                end
            end
        end
    end

    if test_output_after_restart
        @testset "Outputting after a restart" begin 
            @info "Testing outputting after a restart"
            
            for output_format in ["mat", "jld2"]

                folder = "outputs/"
                isdir(folder) && rm(folder, force = true, recursive = true) #make a new folder
                mkdir(folder)

                #run simulation first time with no zipping
                output_test(; dt  = 0.5, 
                            end_time = 10., 
                            output_freq = 0.5, 
                            output_format = output_format, 
                            output_path = folder,
                            zip_format = "nc", 
                            dump_vel = true, 
                            pchkpt_freq = 1.)

                #get the time that the first file was outputted (test for https://github.com/WAVI-ice-sheet-model/WAVI.jl/issues/35)
                first_file_name = joinpath("outputs",string("outfile0000000001.", output_format));
                @test isfile(first_file_name) #check the first output point exists
                dt1 = Dates.unix2datetime(mtime(first_file_name)) #time of last modification

                #run again with different niter0
                sim = output_test(; dt  = 0.5, 
                            end_time = 20., 
                            output_freq = 0.5, 
                            niter0 = 20,
                            output_format = output_format, 
                            output_path = folder,
                            zip_format = "nc", 
                            dump_vel = true, 
                            pchkpt_freq = 1.)

                #check that the time that first file modified has not changed -- i.e. simulation has not touched outputs at earlier times than it should have            
                dt2 = Dates.unix2datetime(mtime(first_file_name)) 
                @test dt1 == dt2
                
                #test variables read from nc file
                fname = string(folder, sim.output_params.prefix, ".nc")
                println(fname)
                NCDataset(fname, "r") do ds
                    x = ds["x"][:]
                    y = ds["y"][:]
                    h = ds["h"][:,:,:]
                    u = ds["u"][:,:,:]
                    v = ds["v"][:,:,:]
                    b = ds["b"][:,:,:]
                    t = ds["TIME"][:]

                    @test ds["TIME"][:] == 0.5:0.5:20.
                    @test size(x) == (80,)
                    @test size(y) == (10,)
                    @test size(b) == (length(x), length(y), length(t))
                    @test size(u) == (length(x), length(y), length(t))
                    @test size(v) == (length(x), length(y), length(t))
                    @test size(h) == (length(x), length(y), length(t))
                end

                #check we have the right number of output files 
                foldersim = sim.output_params.output_path 
                @test foldersim[end] == '/' #test that we do have the / at end of path
                files = [string(foldersim,f) for f in  readdir(foldersim) if endswith(f, output_format)]
                println(files)
                @test length(files) == 40    #check there are the correct number of output files
                
                #now repeat for just one timestep
                sim = output_test(; dt  = 0.5, 
                                end_time = 20.5, 
                                output_freq = 0.5, 
                                niter0 = 40,
                                output_format = output_format, 
                                output_path = folder,
                                zip_format = "nc", 
                                dump_vel = true, 
                                pchkpt_freq = 1.)
                
                #test variables read from nc file
                fname = string(folder, sim.output_params.prefix, ".nc")
                println(fname)

                NCDataset(fname, "r") do ds
                    x = ds["x"][:]
                    y = ds["y"][:]
                    h = ds["h"][:,:,:]
                    u = ds["u"][:,:,:]
                    v = ds["v"][:,:,:]
                    b = ds["b"][:,:,:]
                    t = ds["TIME"][:]

                    @test ds["TIME"][:] == 0.5:0.5:20.5
                    @test size(x) == (80,)
                    @test size(y) == (10,)
                    @test size(b) == (length(x), length(y), length(t))
                    @test size(u) == (length(x), length(y), length(t))
                    @test size(v) == (length(x), length(y), length(t))
                    @test size(h) == (length(x), length(y), length(t))
                end

                #check we have the right number of files
                foldersim = sim.output_params.output_path 
                @test foldersim[end] == '/' #test that we do have the / at end of path
                files = [string(foldersim,f) for f in  readdir(foldersim) if endswith(f, output_format)]
                println(files)
                @test length(files) == 41    #check there are the correct number of output files

                #check we have dumped the velocity
                @test isfile(string(foldersim, "outfile.nc")) #check the zipped file exists


                # delete everything you just made
                rm(folder, force = true, recursive = true)
                foreach(rm, filter(endswith(".mat"), readdir()))
                foreach(rm, filter(endswith(".jld2"), readdir()))
                foreach(rm, filter(endswith(".bin"), readdir()))
            end
        end
    end

    if test_multigrid_output
        @testset "Multi-grid NetCDF output" begin
            @info "Testing multi-grid NetCDF output..."
            
            folder = "outputs/"
            isdir(folder) && rm(folder, force = true, recursive = true)
            mkdir(folder)
            
            # Create a grid and model to get multi-grid variables
            grid = Grid()
            bed = WAVI.mismip_plus_bed
            solver_params = SolverParams(maxiter_picard = 1)
            model = Model(grid = grid, bed_elevation = bed, solver_params = solver_params)
            
            # Create multi-grid variables with different dimensions
            nx, ny = grid.nx, grid.ny
            mg_outputs = (
                h = model.fields.gh.h,                    # h-grid (80×10)
                b = model.fields.gh.b,                    # h-grid (80×10)
                u = model.fields.gh.u,                    # h-grid in this case (80×10)
                v = model.fields.gh.v,                    # h-grid in this case (80×10)
                u_vel = model.fields.gu.u,                # u-grid (81×10) - u-velocity on u-grid
                v_vel = model.fields.gv.v                 # v-grid (80×11) - v-velocity on v-grid
            )
            
            # Run simulation with multi-grid outputs
            timestepping_params = TimesteppingParams(dt = 0.5, end_time = 2.0)
            output_params = OutputParams(
                outputs = mg_outputs,
                output_freq = 0.5,
                output_format = "jld2",
                output_path = folder,
                zip_format = "nc",
                prefix = "multigrid_test",
                output_start = true
            )
            
            simulation = Simulation(
                model = model,
                timestepping_params = timestepping_params,
                output_params = output_params
            )
            
            run_simulation!(simulation)
            
            # Test the resulting NetCDF file
            nc_file = joinpath(folder, "multigrid_test.nc")
            @test isfile(nc_file)
            
            NCDataset(nc_file, "r") do ds
                # Test coordinate systems for each grid
                @test haskey(ds, "x")     # h-grid coordinates
                @test haskey(ds, "y")
                @test haskey(ds, "x_u")   # u-grid coordinates
                @test haskey(ds, "y_u")
                @test haskey(ds, "x_v")   # v-grid coordinates
                @test haskey(ds, "y_v")
                
                # Test coordinate dimensions (default 80×10 grid)
                @test length(ds["x"]) == 80      # nx
                @test length(ds["y"]) == 10      # ny
                @test length(ds["x_u"]) == 81    # nx+1
                @test length(ds["y_u"]) == 10    # ny
                @test length(ds["x_v"]) == 80    # nx
                @test length(ds["y_v"]) == 11    # ny+1
                
                # Test variable dimensions match expected grids
                # h-grid variables (80 × 10)
                @test size(ds["h"]) == (80, 10, 5)  # nx × ny × time
                @test size(ds["b"]) == (80, 10, 5)
                @test size(ds["u"]) == (80, 10, 5)  # model u,v are on h-grid
                @test size(ds["v"]) == (80, 10, 5)
                
                # u-grid variables (81 × 10)
                @test size(ds["u_vel"]) == (81, 10, 5)  # (nx+1) × ny × time
                
                # v-grid variables (80 × 11)
                @test size(ds["v_vel"]) == (80, 11, 5)  # nx × (ny+1) × time
                
                # Test time coordinate
                @test ds["TIME"][:] ≈ [0.0, 0.5, 1.0, 1.5, 2.0]
            end
            
            # Clean up
            rm(folder, force = true, recursive = true)
        end
        
        @testset "Grid detection and error handling" begin
            @info "Testing grid detection and error handling..."
            
            folder = "outputs/"
            isdir(folder) && rm(folder, force = true, recursive = true)
            mkdir(folder)
            
            # Create a standard model
            grid = Grid()
            bed = WAVI.mismip_plus_bed
            solver_params = SolverParams(maxiter_picard = 1)
            model = Model(grid = grid, bed_elevation = bed, solver_params = solver_params)
            
            # Create outputs with one invalid dimension variable
            invalid_outputs = (
                h = model.fields.gh.h,      # Valid h-grid (80×10)
                invalid = zeros(12, 9)       # Invalid dimensions
            )
            
            timestepping_params = TimesteppingParams(dt = 1.0, end_time = 1.0)
            output_params = OutputParams(
                outputs = invalid_outputs,
                output_freq = 1.0,
                output_format = "jld2",
                output_path = folder,
                zip_format = "nc",
                prefix = "error_test"
            )
            
            simulation = Simulation(
                model = model,
                timestepping_params = timestepping_params,
                output_params = output_params
            )
            
            # This should warn about invalid grid dimensions but still complete
            @test_logs (:warn, r"found an output variable.*with spatial dimensions.*that do not match any known grid type") run_simulation!(simulation)
            
            # Verify that the NetCDF file was created but only contains valid variables
            nc_file = joinpath(folder, "error_test.nc")
            @test isfile(nc_file)
            
            NCDataset(nc_file, "r") do ds
                # Should contain the valid h variable
                @test haskey(ds, "h")
                # Should NOT contain the invalid variable
                @test !haskey(ds, "invalid")
            end
            
            # Clean up
            rm(folder, force = true, recursive = true)
        end
        
        @testset "Backward compatibility" begin
            @info "Testing backward compatibility with existing h-grid only outputs..."
            
            folder = "outputs/"
            isdir(folder) && rm(folder, force = true, recursive = true)
            mkdir(folder)
            
            # Test using the standard output_test function with only h-grid variables
            sim = output_test(
                dt = 1.0,
                end_time = 2.0,
                output_freq = 1.0,
                output_format = "jld2",
                output_path = folder,
                zip_format = "nc",
                prefix = "compat_test",
                output_start = true
            )
            
            # Test the NetCDF file maintains backward compatibility
            nc_file = joinpath(folder, "compat_test.nc")
            @test isfile(nc_file)
            
            NCDataset(nc_file, "r") do ds
                # Should still have traditional x, y coordinates
                @test haskey(ds, "x")
                @test haskey(ds, "y")
                
                # Should NOT have u/v grid coordinates when only h-grid used
                @test !haskey(ds, "x_u")
                @test !haskey(ds, "y_u")
                @test !haskey(ds, "x_v")
                @test !haskey(ds, "y_v")
                
                # Variables should be on h-grid dimensions (default 80×10)
                @test size(ds["h"]) == (80, 10, 3)  # nx × ny × time (0, 1, 2)
                @test size(ds["b"]) == (80, 10, 3)
                @test size(ds["u"]) == (80, 10, 3)  # model u,v are on h-grid
                @test size(ds["v"]) == (80, 10, 3)
                
                # Coordinate dimensions  
                @test length(ds["x"]) == 80
                @test length(ds["y"]) == 10
            end
            
            # Clean up
            rm(folder, force = true, recursive = true)
        end
    end

    if test_output_errors
        @testset "Outputting errors" begin 
            @info "Testing outputting errors"
            #check error for weird output format 
            @test_throws ArgumentError output_test(output_format = "incorrect_format")

            #check that non-standard zip_format reverts to none
            sim = output_test(end_time = 1., zip_format = "incorrect_zip_format")
            @test sim.output_params.zip_format == "none"

            @test_throws ArgumentError output_test(output_freq =  -2.)

        end
    end
end

