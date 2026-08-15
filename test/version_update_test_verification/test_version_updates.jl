using WAVI, Test, JLD2

@testset "Version Updates" begin
    @info "testing the current version of WAVI against previous version...."

    function version_update_test()
        #Grid and boundary conditions
        nx = 80
        ny = 10
        nσ = 4
        x0 = 0.0
        y0 = -40000.0
        dx = 8000.0
        dy = 8000.0
        h_mask=trues(nx,ny)
        u_iszero = falses(nx+1,ny); u_iszero[1,:].=true
        v_iszero=falses(nx,ny+1); v_iszero[:,1].=true; v_iszero[:,end].=true

        #alternative bc notations
        u_iszero = ["north"]
        v_iszero = ["west", "east"]
        grid = Grid(nx = nx, 
                    ny = ny,   
                    nσ = nσ, 
                    x0 = x0, 
                    y0 = y0, 
                    dx = dx, 
                    dy = dy,
                    h_mask = h_mask, 
                    u_iszero = u_iszero, 
                    v_iszero = v_iszero)

        #Bed 
        bed = WAVI.mismip_plus_bed #function definition

        #solver parameters
        maxiter_picard = 1
        solver_params = SolverParams(maxiter_picard = maxiter_picard)

        #Physical parameters
        default_thickness = 100.0 #set the initial condition this way
        accumulation_rate = 0.3
        default_temperature=265.700709
        params = Params(default_thickness = default_thickness, 
                        accumulation_rate = accumulation_rate,
                        default_temperature = default_temperature)

        #make the model
        model = Model(grid = grid,
                        bed_elevation = bed, 
                        params = params, 
                        solver_params = solver_params)

        #timestepping parameters
        niter0 = 0
        dt = 0.1
        end_time = 100.
        timestepping_params = TimesteppingParams(niter0 = niter0, 
                                                dt = dt, 
                                                end_time = end_time)

        #nb no output parameters

        simulation = Simulation(model = model, 
                            timestepping_params = timestepping_params)
                
        #perform the simulation
        run_simulation!(simulation)
        return simulation
    end

    simulation = version_update_test();

    if  VERSION == v"1.8.3"
        filename = joinpath(dirname(@__FILE__), "v1_8_3_MISMIP_100yr_output_8kmres_maxiter1_timesteppt1.jld2")
        example_output = load(filename)
    elseif  VERSION == v"1.6.2"
        filename = joinpath(dirname(@__FILE__), "v1_6_2_MISMIP_100yr_output_8kmres_maxiter1_timesteppt1.jld2")
        example_output = load(filename)
    elseif VERSION == v"1.6.1"
        filename = joinpath(dirname(@__FILE__), "v1_6_1_MISMIP_100yr_output_8kmres_maxiter1_timesteppt1.jld2")
        example_output = load(filename)
    elseif  VERSION == v"1.5.1"
        filename = joinpath(dirname(@__FILE__), "v1_5_1_MISMIP_100yr_output_8kmres_maxiter1_timesteppt1.jld2")
        example_output = load(filename)
    else
        filename = joinpath(dirname(@__FILE__), "v1_8_3_MISMIP_100yr_output_8kmres_maxiter1_timesteppt1.jld2")
        example_output = load(filename)
    end

    @testset "Approximate comparison" begin
        @test simulation.model.fields.gh.h ≈ example_output["h"]
        @test_broken simulation.model.fields.gu.u ≈ example_output["u"]
        @test_broken simulation.model.fields.gv.v ≈ example_output["v"]
        @test_broken simulation.model.fields.gh.ηav ≈ example_output["viscosity"]
        @test_broken simulation.model.fields.gh.grounded_fraction ≈ example_output["grounded_fraction"]
        @test_broken simulation.model.fields.gh.bed_speed ≈ example_output["bed_speed"]
        rtoltest = 1e-5
        @test rtoltest > maximum(abs.(simulation.model.fields.gh.h .- example_output["h"]))./maximum(abs.(simulation.model.fields.gh.h))
        @test rtoltest > maximum(abs.(simulation.model.fields.gu.u .- example_output["u"]))./maximum(abs.(simulation.model.fields.gu.u))
        @test rtoltest > maximum(abs.(simulation.model.fields.gv.v .- example_output["v"]))./maximum(abs.(simulation.model.fields.gv.v))
        @test rtoltest > maximum(abs.(simulation.model.fields.gh.ηav .- example_output["viscosity"]))./maximum(abs.(simulation.model.fields.gh.ηav))
        @test rtoltest > maximum(abs.(simulation.model.fields.gh.grounded_fraction .- example_output["grounded_fraction"]))./maximum(abs.(simulation.model.fields.gh.grounded_fraction))
        @test rtoltest > maximum(abs.(simulation.model.fields.gh.bed_speed .- example_output["bed_speed"]))./maximum(abs.(simulation.model.fields.gh.bed_speed))
    end

    @testset "Exact comparison" begin 
        # If floating point rounding has changed since reference runs these will be broken.
        @test_broken simulation.model.fields.gh.h == example_output["h"]
        @test_broken simulation.model.fields.gu.u == example_output["u"]
        @test_broken simulation.model.fields.gv.v == example_output["v"]
        @test_broken simulation.model.fields.gh.ηav == example_output["viscosity"]
        @test_broken simulation.model.fields.gh.grounded_fraction == example_output["grounded_fraction"]
        @test_broken simulation.model.fields.gh.bed_speed == example_output["bed_speed"]
    end

end


