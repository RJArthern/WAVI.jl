using Test, WAVI

@testset "MISMIP+ verification experiments" begin 
    @testset "MISMIP+ Ice0 verification experiments" begin
        @info "Performing MISMIP+ Ice0 verification experiments"
        include("MISMIP_PLUS_Ice0.jl")
        simulation=MISMIP_PLUS_test()
        glx=WAVI.get_glx(simulation.model)
        glxtest=glx[[1,div(simulation.model.grid.ny,2),div(simulation.model.grid.ny,2)+1,simulation.model.grid.ny]]
        @test length(glx) == simluation.model.grid.ny #check that the grounding line covers the whole domain in the y-direction
        @test (glxtest[4]-glxtest[1])/(glxtest[4]+glxtest[1]) < 1e-4
        @test (glxtest[2]-glxtest[3])/(glxtest[2]+glxtest[3]) < 1e-4
        @test 480000<glxtest[1]<540000
        @test 480000<glxtest[4]<540000
        @test 430000<glxtest[2]<460000
        @test 430000<glxtest[3]<460000
    end
end