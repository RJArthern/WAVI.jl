using Test, WAVI

@testset "Models" begin
    @info "Testing models...."

    @testset "testing model construction" begin 
    @info "Testing generic model constuction..."
    grid = Grid()
    bed_elevation = zeros(grid.nx, grid.ny)
    model = Model(grid = grid, bed_elevation = bed_elevation)
    @test model isa Model

    #check that a scalar weertman c passed to model emerges as an array
    model = Model(grid = grid, bed_elevation = bed_elevation, params = Params(weertman_c = 1.0))
    @test model.params.weertman_c == 1.0 *ones(grid.nx, grid.ny)

    #check that an array weertman_c works
    model = Model(grid = grid, bed_elevation = bed_elevation, params = Params(weertman_c =1.0*ones(grid.nx, grid.ny)))
    @test model isa Model

    #check that a scalar weertman c passed to model emerges as an array
    model = Model(grid = grid, bed_elevation = bed_elevation, params = Params(accumulation_rate = 1.0))
    @test model.params.accumulation_rate == 1.0 *ones(grid.nx, grid.ny)

    #check that an array weertman_c works
    model = Model(grid = grid, bed_elevation = bed_elevation, params = Params(accumulation_rate =1.0*ones(grid.nx, grid.ny)))
    @test model isa Model
    end

    @testset "testing model construction errors"  begin
    @info "Testing model construction errors"
    #test dimnesion mismatch if size of input weertman_c incompatible
    grid = Grid()
    bed_elevation = zeros(grid.nx, grid.ny)
    @test_throws DimensionMismatch Model(grid = grid, bed_elevation = bed_elevation, params = Params(weertman_c =1.0*ones(grid.nx -1, grid.ny -1)))

    end

    @testset "Testing generic velocity solve on a model" begin
    @info "Testing generic model velocity solve..."
    grid = Grid()
    bed_elevation = zeros(grid.nx, grid.ny)
    model = Model(grid = grid, bed_elevation = bed_elevation)
    update_state!(model)
    @test model isa Model
    end

end
