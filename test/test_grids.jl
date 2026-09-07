using Test, WAVI

@testset "Grids" begin
    @info "Testing Grids...."

    @testset "Testing Grid construction" begin 
    @info "Testing Grid construction..."

    #no arguments
    nx = 80
    ny = 10
    grid = Grid(nx = nx, ny = ny)
    @test grid isa Grid

    #test that nothing default arguments picked up
    @test grid.h_mask == trues(nx,ny)
    @test grid.u_iszero == falses(nx+1,ny)
    @test grid.v_iszero == falses(nx,ny+1)
    @test grid.quadrature_weights == 0.5*[ grid.σ[2] .- grid.σ[1] ; grid.σ[3:end] .- grid.σ[1:end-2] ; grid.σ[end] .- grid.σ[end-1] ] 
    @test grid.quadrature_weights ≈ [0.5;ones(grid.nσ-2);0.5]/(grid.nσ-1)

    #test the string boundary conditions method
    grid = Grid(nx = nx, ny = ny, u_iszero = ["north"])
    @test all(grid.u_iszero[1,:] .== 1)
    grid = Grid(nx = nx, ny = ny, u_iszero = ["south"])
    @test all(grid.u_iszero[end,:] .== 1)
    grid = Grid(nx = nx, ny = ny, u_iszero = ["west"])
    @test all(grid.u_iszero[:,1] .== 1)
    grid = Grid(nx = nx, ny = ny, u_iszero = ["east"])
    @test all(grid.u_iszero[:,end] .== 1)
    grid = Grid(nx = nx, ny = ny, v_iszero = ["North"]) #with upper case
    @test all(grid.v_iszero[1,:] .== 1)
    grid = Grid(nx = nx, ny = ny, v_iszero = ["SouTh"])
    @test all(grid.v_iszero[end,:] .== 1)
    grid = Grid(nx = nx, ny = ny, v_iszero = ["WeSt"])
    @test all(grid.v_iszero[:,1] .== 1)
    grid = Grid(nx = nx, ny = ny, v_iszero = ["East"])
    @test all(grid.v_iszero[:,end] .== 1)
    end

    @testset "Testing Grid construction" begin 
        @info "Testing Grid construction errors..."
    
        #h_mask size incorrect
        @test_throws DimensionMismatch Grid(nx = 10, ny = 10, h_mask = trues(5,5))
        @test_throws DimensionMismatch Grid(nx = 5, ny = 10, h_mask = trues(5,5))
        @test_throws DimensionMismatch Grid(nx = 10, ny = 5, h_mask = trues(5,5))
        nσ = 10
        @test_throws DimensionMismatch Grid(nσ = nσ, quadrature_weights =  [0.5;ones(nσ-1);0.5]/(nσ-1)) #quadrature weights one short
        
        #u_iszero size incorrect
        @test_throws DimensionMismatch Grid(nx = 10, ny = 10, u_iszero = ones(10,10))

        #v_iszero size incorrect
        @test_throws DimensionMismatch Grid(nx = 10, ny = 10, v_iszero = ones(10,10))

        #h_mask non-boolean
        @test_throws ArgumentError Grid(nx = 10, ny = 10, h_mask = 2.0*ones(10,10))

        #u_iszero non-boolean
        @test_throws ArgumentError Grid(nx = 10, ny = 10, u_iszero = 2.0*ones(11,10))

        #v_iszero non-boolean
        @test_throws ArgumentError Grid(nx = 10, ny = 10, v_iszero = 2.0*ones(10,11))

        #non-positive integer number of grid points
        @test_throws ArgumentError Grid(nx = 10.,)
        @test_throws ArgumentError Grid(ny = 10.,)
        @test_throws ArgumentError Grid(nσ = 10.,)
        @test_throws ArgumentError Grid(nx = -5)
        @test_throws ArgumentError Grid(ny = -5)
        @test_throws ArgumentError Grid(nσ = -5)

    end

    @testset "spatial_on_subdomain" begin
        grid = Grid(nx = 8, ny = 6)
        bounds = (2, 5, 3, 6)
        @test WAVI.Grids.spatial_on_subdomain(0.3, grid, bounds) == 0.3
        @test WAVI.Grids.spatial_on_subdomain(nothing, grid, bounds) === nothing
        f = (x, y) -> x + y
        @test WAVI.Grids.spatial_on_subdomain(f, grid, bounds) === f
        a = reshape(collect(1.0:48.0), 8, 6)
        sliced = WAVI.Grids.spatial_on_subdomain(a, grid, bounds)
        @test sliced == a[2:5, 3:6]
        sliced[1, 1] = -1.0
        @test a[2, 3] != -1.0
        a3 = reshape(collect(1.0:96.0), 8, 6, 2)
        sliced3 = WAVI.Grids.spatial_on_subdomain(a3, grid, bounds)
        @test sliced3 == a3[2:5, 3:6, :]
        sliced3[1, 1, 1] = -1.0
        @test a3[2, 3, 1] != -1.0
    end
end
