using Test
using WAVI
using WAVI.Architectures: device
using KernelAbstractions: KernelAbstractions as KA

@testset "Architectures" begin
    @testset "CPU round-trip" begin
        arch = CPU()
        @test arch isa AbstractArchitecture
        @test array_type(arch) === Array
        @test child_architecture(arch) === arch
        @test device(arch) isa KA.CPU
        @test summary(arch) == "CPU"

        x = rand(2, 3)
        @test on_architecture(arch, x) === x
        @test on_architecture(arch, 1.5) === 1.5
        @test on_architecture(arch, nothing) === nothing
        synchronise(arch)
    end

    cuda_ext = Base.get_extension(WAVI, :WAVICUDAExt)
    if cuda_ext === nothing
        @testset "GPU requires CUDA extension" begin
            @test_throws MethodError GPU()
            @test isempty(WAVI.Architectures.LOADED_GPU_BACKENDS)
            arch = @test_logs (:warn, r"No GPU backend loaded") gpu_device()
            @test arch isa CPU
            @test_throws ErrorException gpu_device(; force=true)
        end
    else
        @testset "GPU constructor with CUDA extension" begin
            gpu = GPU()
            @test gpu isa GPU
            @test array_type(gpu) !== Array
            @test child_architecture(gpu) === gpu
            @test summary(gpu) isa String

            detected = gpu_device()
            @test detected isa Union{CPU, GPU}
            if detected isa GPU
                x = [1.0 2.0; 3.0 4.0]
                y = on_architecture(detected, x)
                @test on_architecture(CPU(), y) == x
                synchronise(detected)
            end
        end
    end
end
