using Test, WAVI

@testset "Melting" begin
    @info "Testing melt rate phsyics...."

    @testset "test analytic melt rate construction" begin 
        @info "Testing analytic melt rate construction"
        grid = Grid()
        bed_elevation = zeros(grid.nx, grid.ny)
        model = Model(grid = grid, bed_elevation = bed_elevation)

        #melt rate model
        m1(h) = h
        arguments = (h = model.fields.gh.h,)
        add_melt_rate_model!(model, AnalyticMeltRate(melt_rate_function = m1 , function_arguments = arguments))
        @test model.extra_physics["melt_rate_model"] isa WAVI.AbstractMeltRateModel

        #check that we get a warning and an overwrite if we add another melt rate model
        m2(x) = 2 .* x
        @test_logs (:info, "Model already contained a melt rate model. Overwritten to that just specified...") add_melt_rate_model!(model, AnalyticMeltRate(melt_rate_function = m2,function_arguments = arguments))
        @test model.extra_physics["melt_rate_model"].melt_rate_function(1) == 2
        
    end
    @testset "test analytic melt rate construction errors" begin
        @info "Testing analytic melt rate construction errors"
        grid = Grid()
        bed_elevation = zeros(grid.nx, grid.ny)
        model = Model(grid = grid, bed_elevation = bed_elevation)
        m1(h) = h

        @test_throws ArgumentError AnalyticMeltRate()
        @test_throws ArgumentError AnalyticMeltRate(function_arguments = (h = model.fields.gh.h))
        @test_throws ArgumentError AnalyticMeltRate(melt_rate_function = m1)
        @test_throws ArgumentError AnalyticMeltRate(melt_rate_function = m1, function_arguments = (h = model.fields.gh.h, x = 1) )

        #check that melt rate model not added when grid dimensions incompatible
        m_err(x,y) = x .+ y
        arguments = (draft = zeros(5,5), cavity_thickness = zeros(5,5))
        model = Model(grid = grid, bed_elevation = bed_elevation);
        add_melt_rate_model!(model, AnalyticMeltRate(melt_rate_function = m_err,function_arguments = arguments))
        @test_throws KeyError model.extra_physics["melt_rate_model"]
    end

    @testset "test updates of analytic melt rate" begin 
        @info "Testing analytic melt rate update"
        grid = Grid()
        bed_elevation = -900 .* ones(grid.nx, grid.ny) #low bed so we can make it float everywhere
        initial_conditions = InitialConditions(initial_thickness = 100.0*ones(grid.nx, grid.ny)) #initialize to 100m thickn
        model = Model(grid = grid, bed_elevation = bed_elevation, initial_conditions = initial_conditions)
        m(x) = x
        melt_rate_model = AnalyticMeltRate(melt_rate_function = m,function_arguments = (x = model.fields.gh.h,)) #melt model sending the melt rate to the thickness
        @test all(melt_rate_model.melt_rate .== 0) #any melt rate model should initialze with zero melt
        
        add_melt_rate_model!(model, melt_rate_model) 
        @test all(model.extra_physics["melt_rate_model"].melt_rate .== 100) #adding melt rate model to the model will update the melt rate

        model.fields.gh.h .= 400 #change the thickness artifically
        WAVI.update_melt_rate_model!(model.extra_physics["melt_rate_model"],model)
        @test all(model.extra_physics["melt_rate_model"].melt_rate .== 400) 
        
        #check that this then passes through to the model gh field
        WAVI.update_basal_melt!(model)
        @test all(model.fields.gh.basal_melt .== 400)

        #check that update_basal_melt! can do this on it's own 
        model.fields.gh.h .= 200 #change the thickness artifically
        WAVI.update_basal_melt!(model)
        @test all(model.fields.gh.basal_melt .== 200) 


    end

    @testset "test binary input file melt rate construction" begin 

        @info "Testing binary input file construction and update"
        grid = Grid()
        bed_elevation = -900 .* ones(grid.nx, grid.ny) #low bed so we can make it float everywhere
        initial_conditions = InitialConditions(initial_thickness = 100.0*ones(grid.nx, grid.ny)) #initialize to 100m thickn
        model = Model(grid = grid, bed_elevation = bed_elevation, initial_conditions = initial_conditions)

        #write a binary file 
        m = ones(grid.nx, grid.ny)
        m .= hton.(m)
        filename =  "melt_test_file.bin"
        mfileID =  open(filename,"w")
          write(mfileID, m[:,:])
        close(mfileID) 

        #create a melt rate model from this input
        binfile_melt_rate = BinfileMeltRate(input_filename =filename,nx = grid.nx, ny = grid.ny)
        @test binfile_melt_rate isa WAVI.BinfileMeltRate

        #add to the model
        add_melt_rate_model!(model, binfile_melt_rate) 
        @test model.extra_physics["melt_rate_model"] isa WAVI.BinfileMeltRate
        @test all(model.extra_physics["melt_rate_model"].melt_rate .== 1)
        @test all(binfile_melt_rate.melt_rate .== 1)

        #change the binary file and check we can update
        m = 2 .* ones(grid.nx, grid.ny)
        m .= hton.(m)
        mfileID =  open(filename,"w")
          write(mfileID, m[:,:])
        close(mfileID)
        update_state!(model)
        @test all(model.extra_physics["melt_rate_model"].melt_rate .== 2)
        @test all(binfile_melt_rate.melt_rate .== 2)

        rm(filename)

        
    end

    @testset "test binary input file melt rate construction errors" begin 

        @info "Testing input file melt rate construction errors"
        @test_throws ArgumentError BinfileMeltRate()
        @test_throws ArgumentError BinfileMeltRate(input_filename = "a_file_that_does_not_exist.bin")

        #size incompatibilities
        grid = Grid()
        bed_elevation = ones(grid.nx, grid.ny) #low bed so we can make it float everywhere
        model = Model(grid = grid, bed_elevation = bed_elevation)
        m = ones(grid.nx+1, grid.ny+1)
        m .= hton.(m)
        filename =  "melt_test_file.bin"
        mfileID =  open(filename,"w")
          write(mfileID, m[:,:])
        close(mfileID) 
        #check that we can't add file of the wrong size to a BinfileMeltRate
        @test_throws DimensionMismatch BinfileMeltRate(input_filename = filename, nx = grid.nx, ny = grid.ny)
        
        #check that we can't add a melt rate model of the wrong size to model
        binfile_melt_rate = BinfileMeltRate(input_filename = filename, nx = grid.nx+1, ny = grid.ny+1)
        add_melt_rate_model!(model, binfile_melt_rate)
        @test_throws KeyError model.extra_physics["melt_rate_model"]

        rm(filename)
        
    end


end
