using Test, WAVI

@testset "Simulation" begin
    @info "Testing simulations..."
    @testset "Simulation construction" begin 
    @info "testing simulation construction"

    grid = Grid()
    bed_elevation = zeros(grid.nx, grid.ny)
    model = Model(grid = grid, bed_elevation = bed_elevation)
    simulation = Simulation(model = model, timestepping_params = TimesteppingParams(end_time = 1.))
    @test simulation isa Simulation
    end

    @testset "Timestepping" begin 
    @info "testing generic timestep"

    grid = Grid()
    bed_elevation = zeros(grid.nx, grid.ny)
    model = Model(grid = grid, bed_elevation = bed_elevation)
    simulation = Simulation(model = model, timestepping_params = TimesteppingParams(end_time = 1.))
    timestep!(simulation)
    @test simulation isa Simulation
    end


end
