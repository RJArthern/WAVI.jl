using Test, WAVI, LinearAlgebra

@testset "Fields" begin
    @info "Testing Fields...."

    @testset "Testing HGrid" begin 
    @info "Testing HGrid construction..."
    hgrid = WAVI.HGrid(nxh = 10, nyh = 10, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(10,10))
    @test hgrid isa WAVI.HGrid
    end

    @testset "Testing HGrid errors" begin 
        @info "Testing HGrid size input errors..."
        @test_throws DimensionMismatch WAVI.HGrid(nxh = 11, nyh = 10, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(10,10))
        @test_throws DimensionMismatch WAVI.HGrid(nxh = 10, nyh = 11, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(10,10))
        @test_throws DimensionMismatch WAVI.HGrid(nxh = 10, nyh = 10, mask = trues(11,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(10,10))
        @test_throws DimensionMismatch WAVI.HGrid(nxh = 10, nyh = 10, mask = trues(10,10), h_isfixed = falses(11,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(10,10))
        @test_throws DimensionMismatch WAVI.HGrid(nxh = 10, nyh = 10, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(11,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(10,10))
        @test_throws DimensionMismatch WAVI.HGrid(nxh = 10, nyh = 10, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(11,10), ηav = ones(10,10), grounded_fraction = ones(10,10))
        @test_throws DimensionMismatch WAVI.HGrid(nxh = 10, nyh = 10, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(11,10), grounded_fraction = ones(10,10))
        @test_throws DimensionMismatch WAVI.HGrid(nxh = 10, nyh = 10, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(11,10))


    end

    @testset "Testing UGrid" begin 
        @info "Testing UGrid construction..."
        ugrid = WAVI.UGrid(nxu = 10, nyu = 10, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test ugrid isa WAVI.UGrid
    end

    @testset "Testing UGrid errors" begin 
        @info "Testing UGrid size input errors..."
        @test_throws DimensionMismatch WAVI.UGrid(nxu = 11, nyu = 10, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test_throws DimensionMismatch WAVI.UGrid(nxu = 10, nyu = 11, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test_throws DimensionMismatch WAVI.UGrid(nxu = 10, nyu = 10, mask = trues(11,10), levels = 5, dx = 10., dy = 10.)
    end
    
    @testset "Testing VGrid" begin 
        @info "Testing VGrid construction..."
        vgrid = WAVI.VGrid(nxv = 10, nyv = 10, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test vgrid isa WAVI.VGrid
    end

    @testset "Testing UGrid errors" begin 
        @info "Testing UGrid size input errors..."
        @test_throws DimensionMismatch WAVI.VGrid(nxv = 11, nyv = 10, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test_throws DimensionMismatch WAVI.VGrid(nxv = 10, nyv = 11, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test_throws DimensionMismatch WAVI.VGrid(nxv = 10, nyv = 10, mask = trues(11,10), levels = 5, dx = 10., dy = 10.)
    end

    @testset "Testing CGrid" begin 
        @info "Testing CGrid construction..."
        cgrid = WAVI.CGrid(nxc = 10, nyc = 10, mask = trues(10,10))
        @test cgrid isa WAVI.CGrid
    end

    @testset "Testing CGrid errors" begin 
        @info "Testing CGrid size input errors..."
        @test_throws DimensionMismatch WAVI.CGrid(nxc = 11, nyc = 10, mask = trues(10,10))
        @test_throws DimensionMismatch WAVI.CGrid(nxc = 10, nyc = 11, mask = trues(10,10))
        @test_throws DimensionMismatch WAVI.CGrid(nxc = 10, nyc = 10, mask = trues(11,10))
    end

    @testset "Testing SigmaGrid" begin 
        @info "Testing SigmaGrid construction..."
        params = Params()
        grid = Grid()
        sigmagrid = WAVI.SigmaGrid(nxs = 10, nys = 10, nσs = 10, σ = collect(range(0., 1., length = 10)), η = ones(10,10,10), θ = ones(10,10,10), Φ = ones(10,10,10), glen_b = fill(WAVI.glen_b(params.default_temperature,params.default_damage,params.glen_a_ref,params.glen_n, params.glen_a_activation_energy, params.glen_temperature_ref, params.gas_const),10,10,10))
        @test sigmagrid isa WAVI.SigmaGrid
        @test dot(sigmagrid.quadrature_weights,ones(eltype(sigmagrid.quadrature_weights),size(sigmagrid.quadrature_weights))) ≈ 1.0
        @test dot(sigmagrid.quadrature_weights,sigmagrid.σ) ≈ 0.5
        sigmagridirregular = WAVI.SigmaGrid(nxs = 10, nys = 10, nσs = 10, σ = collect(range(0., 1., length = 10)).^2, η = ones(10,10,10), θ = ones(10,10,10), Φ = ones(10,10,10), glen_b = fill(WAVI.glen_b(params.default_temperature,params.default_damage,params.glen_a_ref,params.glen_n, params.glen_a_activation_energy, params.glen_temperature_ref, params.gas_const),10,10,10))
        @test sigmagridirregular isa WAVI.SigmaGrid
        @test dot(sigmagridirregular.quadrature_weights,ones(eltype(sigmagridirregular.quadrature_weights),size(sigmagridirregular.quadrature_weights))) ≈ 1.0
        @test dot(sigmagridirregular.quadrature_weights,sigmagridirregular.σ) ≈ 0.5
    end

end
