using Test, WAVI

@testset "TimesteppingParams" begin
    @info "Testing TimesteppingParams..."
    @testset "TimesteppingParams construction" begin 
    @info "testing TimesteppingParams construction"
    #passing end time
    timestepping_params = TimesteppingParams(dt = 0.5, end_time = 100.)
    @test timestepping_params isa TimesteppingParams
    @test timestepping_params.n_iter_total == 200
    #passing total iterations
    timestepping_params = TimesteppingParams(dt = 0.5, n_iter_total = 100)
    @test timestepping_params isa TimesteppingParams
    @test timestepping_params.end_time == 50
    #passing both end time and total iterations
    timestepping_params = TimesteppingParams(dt = 0.5, n_iter_total = 100, end_time = 50.) 
    @test timestepping_params isa TimesteppingParams
    end

    @testset "TimesteppingParams errors" begin 
        @info "testing TimesteppingParams errors"
        @test_throws ArgumentError TimesteppingParams(dt = -0.5) #negative timestep
        @test_throws ArgumentError TimesteppingParams(end_time = -1.) #negative end time
        @test_throws ArgumentError TimesteppingParams(n_iter_total = -1) #negative total iterations 
        @test_throws ArgumentError TimesteppingParams(end_time = 1., dt = 1., n_iter_total = 10) #incompatible endtime and total iters
        @test_throws ArgumentError TimesteppingParams(dt = 0.5) #pass neither of end time or total iterations --> throw error
        @test_throws ArgumentError TimesteppingParams(dt = 0.5, end_time = Inf) 
        @test_throws ArgumentError TimesteppingParams(dt = 0.5, ntimesteps_velocity_update = 0.5) #non integer sub sampling step
        
    end
end
