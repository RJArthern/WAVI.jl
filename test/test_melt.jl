using Test, WAVI

@testset "Melting" begin
    @info "Testing melt rate phsyics...."

    @testset "test analytic melt rate" begin 
    @info "Testing analytic melt rate construction"
    grid = Grid()
    bed_elevation = zeros(grid.nx, grid.ny)
    model = Model(grid = grid, bed_elevation = bed_elevation)

    #melt rate model
    m1(draft,cavity_thickness) = 0.2.* tanh.(cavity_thickness/75).* max.((-100 .- draft), 0)
    ρi = 918.0
    ρw = 1028.0
    draft = -(ρi / ρw) * model.fields.gh.h
    cavity_thickness = (draft .- model.fields.gh.b)
    arguments = (draft = draft, cavity_thickness = cavity_thickness)
    add_melt_rate_model!(model, AnalyticMeltRate(melt_rate_function = m1,function_arguments = arguments))
    @test model.extra_physics["melt_rate_model"] isa WAVI.AbstractMeltRateModel

    #check that we get a warning if we try to add another melt rate model
    @test_logs (:info, "Model already contained a melt rate model. Overwritten to that just specified...") add_melt_rate_model!(model, AnalyticMeltRate(melt_rate_function = m1,function_arguments = arguments))


    @info "Testing analytic melt rate construction errors"
    @test_throws ArgumentError AnalyticMeltRate()
    @test_throws ArgumentError AnalyticMeltRate(function_arguments = arguments)
    @test_throws ArgumentError AnalyticMeltRate(melt_rate_function = m1)
    @test_throws ArgumentError AnalyticMeltRate(melt_rate_function = m1, function_arguments = (draft = draft, cavity_thickness = cavity_thickness, x = 1) )

    #check that melt rate model not added when grid dimensions incompatible
    arguments = (draft = zeros(5,5), cavity_thickness = zeros(5,5))
    model = Model(grid = grid, bed_elevation = bed_elevation);
    add_melt_rate_model!(model, AnalyticMeltRate(melt_rate_function = m1,function_arguments = arguments))
    @test_throws KeyError model.extra_physics["melt_rate_model"]
    end

    @testset "test updates of analytic melt rate" begin 
        @info "Testing analytic melt rate update"
        grid = Grid()
        bed_elevation = zeros(grid.nx, grid.ny)
        model = Model(grid = grid, bed_elevation = bed_elevation)
        add_melt_rate_model!(model, AnalyticMeltRate(melt_rate_function = (x)->x,function_arguments = (x = model.fields.gh.h,))) #melt model sending the melt rate to the thickness
        model.fields.gh.h .= 400 #change the thickness artifically
        @test all(model.extra_physics["melt_rate_model"].melt_rate .== 400) #use of pointers means the melt rate should also update
        #...but we can also do it explicity
        WAVI.update_melt_rate_model!(model.extra_physics["melt_rate_model"],model)
        @test all(model.extra_physics["melt_rate_model"].melt_rate .== 400) 
        
        #check that this then passes through to the model gh field
        WAVI.update_basal_melt!(model)
        @test all(model.fields.gh.basal_melt .== 400)


    end
end
