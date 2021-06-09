using Test, WAVI

@testset "Grids" begin
    @info "Testing Grids...."

    @testset "Testing Grid construction" begin 
    @info "Testing Grid constuction..."

    #no arguments
    nx = 80
    ny = 10
    grid = Grid(nx = nx, ny = ny)
    @test grid isa Grid

    #test that nothing default arguments picked up
    @test grid.h_mask == trues(nx,ny)
    @test grid.u_iszero == falses(nx+1,ny)
    @test grid.v_iszero == falses(nx,ny+1)
    @test grid.quadrature_weights ==  [0.5;ones(grid.nσ-2);0.5]/(grid.nσ-1)

    end

    @testset "Testing Grid construction" begin 
        @info "Testing Grid constuction errors..."
    
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
end