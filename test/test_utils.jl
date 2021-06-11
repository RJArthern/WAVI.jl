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
end
