using Test, WAVI

@testset "Melting" begin
    @testset "test binary input file melt rate construction and update" begin 

        @info "Testing binary input file construction and update"
        grid = Grid()
        bed_elevation = -900 .* ones(grid.nx, grid.ny) #low bed so we can make it float everywhere
        initial_conditions = InitialConditions(initial_thickness = 100.0*ones(grid.nx, grid.ny)) #initialize to 100m thickn

        #write a binary file 
        m = ones(grid.nx, grid.ny)
        m .= hton.(m)
        filename =  "melt_test_file.bin"
        mfileID =  open(filename,"w")
          write(mfileID, m[:,:])
        close(mfileID) 

        #create a melt rate model from this input
        binfile_melt_rate = BinfileMeltRate(input_filename =filename)
        @test binfile_melt_rate isa BinfileMeltRate

        #create model 
        model = Model(grid = grid, 
                      bed_elevation = bed_elevation, 
                      initial_conditions = initial_conditions,
                      melt_rate = binfile_melt_rate,
                      solver_params = SolverParams(maxiter_picard = 1))
        @test model.melt_rate isa BinfileMeltRate

        #check that melt rate read OK
        update_state!(model)
        @test all(model.fields.gh.basal_melt .== 1)

        #change the binary file and check we can update
        m = 2 .* ones(grid.nx, grid.ny)
        m .= hton.(m)
        mfileID =  open(filename,"w")
          write(mfileID, m[:,:])
        close(mfileID)
        update_state!(model)
        @test all(model.fields.gh.basal_melt .== 2)

        #change the binary file to one of the wrong size, and check we return an error when updating
        m = 2 .* ones(grid.nx+1, grid.ny+1)
        m .= hton.(m)
        mfileID =  open(filename,"w")
          write(mfileID, m[:,:])
        close(mfileID)
        @test_throws DimensionMismatch update_state!(model)

        #remove the file we just made
        rm(filename)

        
    end


end
