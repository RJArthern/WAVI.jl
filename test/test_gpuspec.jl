using Test
using WAVI
using SparseArrays
using WAVI.Architectures: architecture, zeros_on, array_type
using KernelAbstractions: KernelAbstractions as KA

@testset "GPUSpec" begin
    @testset "CPU helpers" begin
        arch = CPU()
        z = zeros_on(arch, Float64, 2, 3)
        @test z isa Array{Float64,2}
        @test size(z) == (2, 3)
        @test all(iszero, z)
        @test architecture(arch) === arch
        @test architecture(BasicSpec()) isa CPU
    end

    cuda_ext = Base.get_extension(WAVI, :WAVICUDAExt)
    if cuda_ext === nothing
        @testset "GPUSpec needs a GPU backend" begin
            @test_throws Exception GPUSpec()
        end
    else
        @testset "GPUSpec with CUDA extension" begin
            detected = gpu_device()
            if !(detected isa GPU)
                @test_throws Exception GPUSpec()
            else
                spec = GPUSpec()
                @test spec isa GPUSpec
                @test architecture(spec) isa GPU
                @test array_type(architecture(spec)) !== Array

                grid = Grid(nx = 8, ny = 8, dx = 1.0e3, dy = 1.0e3)
                bed = zeros(grid.nx, grid.ny)
                model = Model(grid = grid, bed_elevation = bed, spec = spec)
                @test architecture(model) isa GPU
                AT = array_type(architecture(spec))
                @test model.fields.gh.h isa AT
                @test model.fields.gu.u isa AT
                @test model.fields.wu.wavelets isa AT
                @test model.fields.gh.samp isa SparseMatrixCSC

                s = WAVI.Utilities.stencil_scratch!(model)
                @test KA.get_backend(s.haar_u) === KA.get_backend(model.fields.gu.u)
                @test KA.get_backend(s.rhs) === KA.get_backend(model.fields.gh.h)

                update_state!(model)
                u_gpu = Array(model.fields.gu.u)
                v_gpu = Array(model.fields.gv.v)

                cpu_model = Model(grid = grid, bed_elevation = bed)
                update_state!(cpu_model)
                @test u_gpu ≈ cpu_model.fields.gu.u rtol = 1e-8 atol = 1e-8
                @test v_gpu ≈ cpu_model.fields.gv.v rtol = 1e-8 atol = 1e-8
            end
        end
    end
end
