module WAVICUDAExt

using WAVI
using CUDA: CUDA, CuArray, CUDABackend
import WAVI.Architectures: GPU, CPU, array_type, on_architecture

# No-argument GPU() exists only after `using CUDA`. Core WAVI has no CUDA types.
# Need to think about how this can be done neatly - maybe on top of a SLURM script?
# or, defined in driver file?
# In future, want to be able to construct GPU architectures from other vendor packages.
# E.g. `GPU(AMDGPU.ROCmBackend())`, Intel's oneAPI backend, etc.
GPU() = GPU(CUDABackend(; always_inline=true))

array_type(::GPU{<:CUDABackend}) = CuArray

on_architecture(::GPU{<:CUDABackend}, a::Array) = CuArray(a)
on_architecture(::GPU{<:CUDABackend}, a::CuArray) = a
on_architecture(::CPU, a::CuArray) = Array(a)

function __init__()
    if CUDA.functional()
        WAVI.Architectures.register_gpu_backend!(GPU())
    end
end

end
