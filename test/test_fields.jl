using Test, WAVI, LinearAlgebra, Adapt

using WAVI.Fields: HGrid, UGrid, VGrid, CGrid, SigmaGrid

@testset "Fields" begin
    @info "Testing Fields...."

    @testset "Testing HGrid" begin 
    @info "Testing HGrid construction..."
    hgrid = HGrid(nxh = 10, nyh = 10, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(10,10), preBfactor = ones(10,10))
    @test hgrid isa HGrid
    end

    @testset "Testing HGrid errors" begin 
        @info "Testing HGrid size input errors..."
        @test_throws DimensionMismatch HGrid(nxh = 11, nyh = 10, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(10,10),preBfactor = ones(10,10))
        @test_throws DimensionMismatch HGrid(nxh = 10, nyh = 11, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(10,10),preBfactor = ones(10,10))
        @test_throws DimensionMismatch HGrid(nxh = 10, nyh = 10, mask = trues(11,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(10,10),preBfactor = ones(10,10))
        @test_throws DimensionMismatch HGrid(nxh = 10, nyh = 10, mask = trues(10,10), h_isfixed = falses(11,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(10,10),preBfactor = ones(10,10))
        @test_throws DimensionMismatch HGrid(nxh = 10, nyh = 10, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(11,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(10,10),preBfactor = ones(10,10))
        @test_throws DimensionMismatch HGrid(nxh = 10, nyh = 10, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(11,10), ηav = ones(10,10), grounded_fraction = ones(10,10),preBfactor = ones(10,10))
        @test_throws DimensionMismatch HGrid(nxh = 10, nyh = 10, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(11,10), grounded_fraction = ones(10,10),preBfactor = ones(10,10))
        @test_throws DimensionMismatch HGrid(nxh = 10, nyh = 10, mask = trues(10,10), h_isfixed = falses(10,10), b = ones(10,10), h = ones(10,10), ηav = ones(10,10), grounded_fraction = ones(11,10),preBfactor = ones(10,10))


    end

    @testset "Testing UGrid" begin 
        @info "Testing UGrid construction..."
        ugrid = UGrid(nxu = 10, nyu = 10, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test ugrid isa UGrid
    end

    @testset "Testing UGrid errors" begin 
        @info "Testing UGrid size input errors..."
        @test_throws DimensionMismatch UGrid(nxu = 11, nyu = 10, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test_throws DimensionMismatch UGrid(nxu = 10, nyu = 11, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test_throws DimensionMismatch UGrid(nxu = 10, nyu = 10, mask = trues(11,10), levels = 5, dx = 10., dy = 10.)
    end
    
    @testset "Testing VGrid" begin 
        @info "Testing VGrid construction..."
        vgrid = VGrid(nxv = 10, nyv = 10, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test vgrid isa VGrid
    end

    @testset "Testing UGrid errors" begin 
        @info "Testing UGrid size input errors..."
        @test_throws DimensionMismatch VGrid(nxv = 11, nyv = 10, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test_throws DimensionMismatch VGrid(nxv = 10, nyv = 11, mask = trues(10,10), levels = 5, dx = 10., dy = 10.)
        @test_throws DimensionMismatch VGrid(nxv = 10, nyv = 10, mask = trues(11,10), levels = 5, dx = 10., dy = 10.)
    end

    @testset "Testing CGrid" begin 
        @info "Testing CGrid construction..."
        cgrid = CGrid(nxc = 10, nyc = 10, mask = trues(10,10))
        @test cgrid isa CGrid
    end

    @testset "Testing CGrid errors" begin 
        @info "Testing CGrid size input errors..."
        @test_throws DimensionMismatch CGrid(nxc = 11, nyc = 10, mask = trues(10,10))
        @test_throws DimensionMismatch CGrid(nxc = 10, nyc = 11, mask = trues(10,10))
        @test_throws DimensionMismatch CGrid(nxc = 10, nyc = 10, mask = trues(11,10))
    end

    @testset "Testing SigmaGrid" begin 
        @info "Testing SigmaGrid construction..."
        params = Params()
        grid = Grid()
        sigmagrid = SigmaGrid(nxs = 10, nys = 10, nσs = 10, σ = collect(range(0., 1., length = 10)), η = ones(10,10,10), θ = ones(10,10,10), Φ = ones(10,10,10), strain_history = zeros(10,10,10), glen_b = fill(WAVI.glen_b(params.default_temperature,params.default_damage,params.glen_a_ref,params.glen_n, params.glen_a_activation_energy, params.glen_temperature_ref, params.gas_const),10,10,10))
        @test sigmagrid isa SigmaGrid
        @test dot(sigmagrid.quadrature_weights,ones(eltype(sigmagrid.quadrature_weights),size(sigmagrid.quadrature_weights))) ≈ 1.0
        @test dot(sigmagrid.quadrature_weights,sigmagrid.σ) ≈ 0.5
        sigmagridirregular = SigmaGrid(nxs = 10, nys = 10, nσs = 10, σ = collect(range(0., 1., length = 10)).^2, η = ones(10,10,10), θ = ones(10,10,10), Φ = ones(10,10,10), strain_history = zeros(10,10,10), glen_b = fill(WAVI.glen_b(params.default_temperature,params.default_damage,params.glen_a_ref,params.glen_n, params.glen_a_activation_energy, params.glen_temperature_ref, params.gas_const),10,10,10))
        @test sigmagridirregular isa SigmaGrid
        @test dot(sigmagridirregular.quadrature_weights,ones(eltype(sigmagridirregular.quadrature_weights),size(sigmagridirregular.quadrature_weights))) ≈ 1.0
        @test dot(sigmagridirregular.quadrature_weights,sigmagridirregular.σ) ≈ 0.5
    end

    @testset "Adapt keeps CPU Arrays and host operators" begin
        @info "Testing Adapt of grid fields..."
        hgrid = HGrid(nxh = 8, nyh = 8, mask = trues(8,8), h_isfixed = falses(8,8),
                      b = ones(8,8), h = ones(8,8), ηav = ones(8,8),
                      grounded_fraction = ones(8,8), preBfactor = ones(8,8))
        adapted = Adapt.adapt(Array, hgrid)
        @test adapted.h isa Array
        @test adapted.h == hgrid.h
        @test adapted.samp === hgrid.samp
        @test adapted.crop === hgrid.crop
        @test adapted.cent_xy === hgrid.cent_xy

        ugrid = UGrid(nxu = 8, nyu = 8, mask = trues(8,8), levels = 3, dx = 10.0, dy = 10.0)
        u_adapted = Adapt.adapt(Array, ugrid)
        @test u_adapted.u isa Array
        @test u_adapted.samp === ugrid.samp
    end

end
