using Test, WAVI

@testset "Fields" begin
    @info "Testing Fields...."

    @testset "Testing HGrid" begin 
    @info "Testing HGrid constuction..."
    hgrid = WAVI.HGrid(Nx = 10, Ny = 10, mask = trues(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10))
    @test hgrid isa WAVI.HGrid
    end

    @testset "Testing HGrid errors" begin 
        @info "Testing HGrid size input errors..."
        @test_throws DimensionMismatch WAVI.HGrid(Nx = 11, Ny = 10, mask = trues(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10))
        @test_throws DimensionMismatch WAVI.HGrid(Nx = 10, Ny = 11, mask = trues(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10))
        @test_throws DimensionMismatch WAVI.HGrid(Nx = 10, Ny = 10, mask = trues(11,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10))
        @test_throws DimensionMismatch WAVI.HGrid(Nx = 10, Ny = 10, mask = trues(10,10), b = ones(11,10), h = ones(10,10), ηav = ones(10,10))
        @test_throws DimensionMismatch WAVI.HGrid(Nx = 10, Ny = 10, mask = trues(10,10), b = ones(10,10), h = ones(11,10), ηav = ones(10,10))
        @test_throws DimensionMismatch WAVI.HGrid(Nx = 10, Ny = 10, mask = trues(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(11,10))
    end

    @testset "Testing UGrid" begin 
        @info "Testing UGrid constuction..."
        ugrid = WAVI.UGrid(Nx = 10, Ny = 10, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test ugrid isa WAVI.UGrid
        end

    @testset "Testing UGrid errors" begin 
        @info "Testing UGrid size input errors..."
        @test_throws DimensionMismatch WAVI.UGrid(Nx = 11, Ny = 10, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test_throws DimensionMismatch WAVI.UGrid(Nx = 10, Ny = 11, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test_throws DimensionMismatch WAVI.UGrid(Nx = 10, Ny = 10, mask = trues(11,10), levels = 5, dx = 10., dy = 10.)
    end
    
    @testset "Testing VGrid" begin 
        @info "Testing VGrid constuction..."
        vgrid = WAVI.VGrid(Nx = 10, Ny = 10, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test vgrid isa WAVI.VGrid
        end

    @testset "Testing UGrid errors" begin 
        @info "Testing UGrid size input errors..."
        @test_throws DimensionMismatch WAVI.VGrid(Nx = 11, Ny = 10, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test_throws DimensionMismatch WAVI.VGrid(Nx = 10, Ny = 11, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test_throws DimensionMismatch WAVI.VGrid(Nx = 10, Ny = 10, mask = trues(11,10), levels = 5, dx = 10., dy = 10.)
    end

    @testset "Testing CGrid" begin 
        @info "Testing CGrid constuction..."
        ugrid = WAVI.CGrid(Nx = 10, Ny = 10, mask = trues(10,10))
        @test ugrid isa WAVI.CGrid
        end

    @testset "Testing UGrid errors" begin 
        @info "Testing UGrid size input errors..."
        @test_throws DimensionMismatch WAVI.CGrid(Nx = 11, Ny = 10, mask = trues(10,10))
        @test_throws DimensionMismatch WAVI.CGrid(Nx = 10, Ny = 11, mask = trues(10,10))
        @test_throws DimensionMismatch WAVI.CGrid(Nx = 10, Ny = 10, mask = trues(11,10))
    end
end
