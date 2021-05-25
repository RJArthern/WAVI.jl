using Test, WAVI

@testset "Models" begin
    @info "Testing models...."

    @testset "testing model construction" begin 
    @info "Testing generic model constuction..."
    grid = Grid()
    bed_elevation = zeros(grid.nx, grid.ny)
    model = Model(grid = grid, bed_elevation = bed_elevation)
    @test model isa Model
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
