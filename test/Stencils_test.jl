using WAVI
using Test
using KernelAbstractions

@testset "Stencils" begin
    # Very basic tests to make sure kernels compile and run

    nx, ny = 10, 10
    inp = ones(nx+1, ny+1)
    out = zeros(nx, ny)
    dx_inv = 1.0

    WAVI.Stencils.launch!(WAVI.Stencils._diff_x!, out, inp, dx_inv; ndrange = (nx, ny))
    @test all(out .== 0.0)

    WAVI.Stencils.launch!(WAVI.Stencils._avg_x!, out, inp; ndrange = (nx, ny))
    @test all(out .== 1.0)
end

@testset "Picard stencil applies" begin
    grid = Grid(nx = 8, ny = 6, dx = 1000.0, dy = 1000.0, nσ = 4)
    model = Model(grid = grid, bed_elevation = fill(-500.0, grid.nx, grid.ny))

    WAVI.Processes.update_av_speed!(model)
    WAVI.Processes.update_shelf_strain_rate!(model)
    WAVI.Processes.update_βeff_on_uv_grids!(model)
    WAVI.Processes.update_rheological_operators!(model)
    rhs = WAVI.Processes.get_rhs(model)

    @test all(isfinite, model.fields.gh.av_speed)
    @test all(isfinite, model.fields.gh.shelf_strain_rate)
    @test length(rhs) == model.fields.gu.ni + model.fields.gv.ni
    @test all(isfinite, rhs)
end
