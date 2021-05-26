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
    model.extra_physics["melt_rate_model"] = AnalyticMeltRate(melt_rate_function = m1,
                                                            function_arguments = arguments)
    @test model.extra_physics["melt_rate_model"] isa WAVI.AbstractMeltRateModel

    @info "Testing analytic melt rate construction errors"
    @test_throws ArgumentError AnalyticMeltRate()
    @test_throws ArgumentError AnalyticMeltRate(function_arguments = arguments)
    @test_throws ArgumentError AnalyticMeltRate(melt_rate_function = m1)
    end

end
