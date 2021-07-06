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

    @testset "Initial Conditions" begin 
    @info "Testing initial conditions"
    #passing no initial conditions should revert to default values
    grid = Grid()
    params = Params()
    bed_elevation = zeros(grid.nx, grid.ny)
    model = Model(grid = grid, bed_elevation = bed_elevation)
    @test model.initial_conditions.initial_thickness == params.default_thickness*ones(grid.nx, grid.ny)
    @test model.initial_conditions.initial_damage == params.default_damage*ones(grid.nx, grid.ny, grid.nσ)
    @test model.initial_conditions.initial_viscosity == params.default_viscosity*ones(grid.nx, grid.ny, grid.nσ)
    @test model.initial_conditions.initial_temperature == params.default_temperature*ones(grid.nx, grid.ny, grid.nσ)
    
    #test passing arrays of thickness, damage, temperature, and viscosity
    initial_conditions = InitialConditions(initial_thickness = 3.14159 * ones(grid.nx,grid.ny),
                                            initial_viscosity = 2.7182818 * ones(grid.nx, grid.ny, grid.nσ),
                                            initial_temperature = 1.618034 * ones(grid.nx, grid.ny, grid.nσ),
                                            initial_damage = 1.414213 * ones(grid.nx, grid.ny, grid.nσ))
    model = Model(grid = grid, bed_elevation = bed_elevation, initial_conditions = initial_conditions)
    @test model isa Model
    @test model.initial_conditions.initial_thickness == 3.14159 * ones(grid.nx,grid.ny)
    @test model.initial_conditions.initial_viscosity == 2.7182818 * ones(grid.nx, grid.ny, grid.nσ)
    @test model.initial_conditions.initial_temperature == 1.618034 * ones(grid.nx, grid.ny, grid.nσ)
    @test model.initial_conditions.initial_damage == 1.414213 * ones(grid.nx, grid.ny, grid.nσ)
    
    #check fields successfully passed to fields
    @test model.fields.gh.h == 3.14159 * ones(grid.nx,grid.ny)
    @test model.fields.g3d.η == 2.7182818 * ones(grid.nx, grid.ny, grid.nσ)
    @test model.fields.g3d.θ == 1.618034 * ones(grid.nx, grid.ny, grid.nσ)
    @test model.fields.g3d.Φ ==  1.414213 * ones(grid.nx, grid.ny, grid.nσ)

    end

    @testset "Initial Conditions Errors" begin
    @info "Testing initial conditions errors"
    #size incompatibilities:
    grid = Grid()
    bed_elevation = zeros(grid.nx, grid.ny)
    initial_conditions = InitialConditions(initial_thickness = ones((grid.nx + 1), grid.ny))
    @test_throws DimensionMismatch Model(grid = grid, bed_elevation = bed_elevation, initial_conditions = initial_conditions)

    initial_conditions = InitialConditions(initial_viscosity = ones((grid.nx + 1), grid.ny, grid.nσ))
    @test_throws DimensionMismatch Model(grid = grid, bed_elevation = bed_elevation, initial_conditions = initial_conditions)

    initial_conditions = InitialConditions(initial_temperature = ones((grid.nx + 1), grid.ny, grid.nσ))
    @test_throws DimensionMismatch Model(grid = grid, bed_elevation = bed_elevation, initial_conditions = initial_conditions)

    initial_conditions = InitialConditions(initial_damage = ones((grid.nx + 1), grid.ny, grid.nσ))
    @test_throws DimensionMismatch Model(grid = grid, bed_elevation = bed_elevation, initial_conditions = initial_conditions)

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
