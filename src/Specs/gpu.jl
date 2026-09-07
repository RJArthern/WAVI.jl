import WAVI.Architectures: architecture
using WAVI.Architectures: AbstractArchitecture, CPU, gpu_device

export GPUSpec

"""
    GPUSpec(arch=gpu_device(; force=true))

Run the local model on one GPU. Some aspects still run on CPU for now.

`arch` must be a `GPU`. Load the vendor package first (`using CUDA`).
Pass `GPUSpec` as `model.spec` so spec dispatch is preserved. Do not wrap
a `BasicSpec` with a separate architecture keyword.
"""
struct GPUSpec{A <: AbstractArchitecture} <: AbstractSpec
    arch::A
    function GPUSpec(arch::AbstractArchitecture = gpu_device(; force=true))
        if arch isa CPU
            throw(ArgumentError(
                "GPUSpec needs a GPU architecture. Load a vendor package first, " *
                "for example `using CUDA`.",
            ))
        end
        return new{typeof(arch)}(arch)
    end
end

architecture(spec::GPUSpec) = spec.arch
