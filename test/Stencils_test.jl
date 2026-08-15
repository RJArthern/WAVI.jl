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
