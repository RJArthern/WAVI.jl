using Test, WAVI

@testset "WAVI tests" begin
    @testset "Pos Fraction" begin
        @info "testing pos fraction"
        z=[-1.0 -1.0 -1.0;-1.0 1.0 -1.0;-1.0 -1.0 -1.0]
        pfh,pfu,pfv=WAVI.pos_fraction(z)
        @test pfh == [0.0 0.0 0.0;0.0 0.5 0.0;0.0 0.0 0.0]
        @test pfu == [0.0 0.0 0.0;0.0 0.25 0.0;0.0 0.25 0.0;0.0 0.0 0.0]
        @test pfv == [0.0 0.0 0.0 0.0;0.0 0.25 0.25 0.0;0.0 0.0 0.0 0.0]
    end

    @testset "Pos Fraction Mask" begin
        @info "testing post fraction mask"
        z=[-1.0 -1.0 -1.0;-1.0 1.0 -1.0;-1.0 -1.0 -1.0]
        mask=[false false false; false true false; false false false]
        pfh,pfu,pfv=WAVI.pos_fraction(z,mask=mask)
        @test pfh == [0.0 0.0 0.0;0.0 1.0 0.0;0.0 0.0 0.0]
        @test pfu == [0.0 0.0 0.0;0.0 0.5 0.0;0.0 0.5 0.0;0.0 0.0 0.0]
        @test pfv == [0.0 0.0 0.0 0.0;0.0 0.5 0.5 0.0;0.0 0.0 0.0 0.0]
    end

    @testset "Height and volume above floatation" begin
        @info "testing height and volume above floatation functions"
        grid = Grid(nx = 20, ny = 25, dx = 1000., dy = 1000.)
        bed = -200.0 * ones(grid.nx, grid.ny)
        initial_conditions = InitialConditions(initial_thickness = 500.0*ones(grid.nx, grid.ny))
        model = Model(grid = grid, bed_elevation = bed, initial_conditions = initial_conditions)
        update_state!(model)
        @test volume_above_floatation(model.fields.gh.h, model.fields.gh.b, Ref(model.params), model.grid) ≈ (500.0 - (1028.0/918.0)*(200.0) ) .* 25 .* 20 .* 1000.0 .^2

    end
end
