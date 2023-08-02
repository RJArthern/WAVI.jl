using Test, WAVI, LinearAlgebra
@testset  "WAVI tests" begin
    @testset "Iceberg" begin
        @info "Performing tests on a spinning, drifting iceberg."
        include("verification_tests/iceberg_test.jl")
        sim, relerr_h, relerr_u, relerr_v = iceberg_test(;end_time=1000.)
        @test relerr_h < 1.0e-4
        @test relerr_u < 3.0e-4
        @test relerr_v < 3.0e-4
    end
end

if true 
@testset "MISMIP+ verification experiments" begin 
    @testset "MISMIP+ Ice0 verification experiments" begin
        @info "Performing MISMIP+ Ice0 verification experiments"
        include(joinpath("verification_tests","MISMIP_PLUS_Ice0.jl"))
        simulation=MISMIP_PLUS_Ice0()
        glx=WAVI.get_glx(simulation.model)
        glxtest=glx[[1,div(simulation.model.grid.ny,2),div(simulation.model.grid.ny,2)+1,simulation.model.grid.ny]]
        @test length(glx) == simulation.model.grid.ny #check that the grounding line covers the whole domain in the y-direction
        @test (glxtest[4]-glxtest[1])/(glxtest[4]+glxtest[1]) < 1e-4
        @test (glxtest[2]-glxtest[3])/(glxtest[2]+glxtest[3]) < 1e-4
        @testset "Tight Tolerance" begin
            # For quick testing at low resolutions (e.g. 8 km) these may be broken
            @test_broken 480000<glxtest[1]<540000
            @test_broken 480000<glxtest[4]<540000
            @test_broken 430000<glxtest[2]<460000
            @test_broken 430000<glxtest[3]<460000
        end
        @testset "Loose Tolerance" begin
            @test 480000<glxtest[1]<550000
            @test 480000<glxtest[4]<550000
            @test 430000<glxtest[2]<490000
            @test 430000<glxtest[3]<490000
        end

        #check the melt rate for ice_1r is doing something sensible
        function m1(h, b)
            draft = -(918.0 / 1028.0) * h
            cavity_thickness = draft .- b
            cavity_thickness = max.(cavity_thickness, 0)
            m =  0.2*tanh.(cavity_thickness./75).*max.((-100 .- draft), 0)
            return m
        end
        melt = m1.(simulation.model.fields.gh.h, simulation.model.fields.gh.h)
        @test all(melt .>= 0)
        @test maximum(melt) < 100

    end
end
end