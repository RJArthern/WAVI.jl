using Test, WAVI, MAT, JLD2

function output_test(; dt  = 0.5, end_time = 100., output_freq = 5., output_format = "jld2", zip_format = "none", prefix = "outfile", dump_vel = false, output_path = ".")
    grid = Grid() #default grid with nx = 80, ny = 10
    bed = WAVI.mismip_plus_bed
    solver_params = SolverParams(maxiter_picard = 1)
    model = Model(grid = grid, bed_elevation = bed, solver_params = solver_params)
    timestepping_params = TimesteppingParams(dt = dt, end_time = end_time)
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
                        dump_vel = dump_vel)

    simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params,
                        output_params = output_params)
            
    #perform the simulation
    run_simulation!(simulation)
    return simulation
end

@testset "Outputting" begin
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
    @testset "Zipping output" begin 
        @info "Testing zipping output..."
    
        for output_format in ["mat", "jld2"]

        folder = "outputs/"
        isdir(folder) && rm(folder, force = true, recursive = true) #make a new folder
        mkdir(folder)

        sim = output_test(zip_format = "nc", output_path = folder, output_format = "mat")
        
        @test sim isa Simulation
        fname = string(folder, sim.output_params.prefix, ".nc")
        @test isfile(fname) #check the zipped file exists

        #test variables read from nc file
        x = ncread(fname, "x");
        y = ncread(fname, "y");
        h = ncread(fname, "h");
        u = ncread(fname, "u");
        v = ncread(fname, "v")
        b = ncread(fname, "b")
        @test ncread(fname, "TIME") == 5.:5.:100.
        @test size(x) == (80,)
        @test size(y) == (10,)
        @test size(b) == (length(x),length(y),length(t))
        @test size(u) == (length(x),length(y),length(t))
        @test size(v) == (length(x),length(y),length(t))
        @test size(h) == (length(x),length(y),length(t))
 
        #delete the folder
        rm(folder, force = true, recursive = true)
        end

    end
end