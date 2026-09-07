module WAVICUDAExt

using WAVI
using CUDA: CUDA, CuArray, CUDABackend
import WAVI.Architectures: GPU, CPU, array_type, on_architecture, assign_local_device!

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

"""
    assign_local_device!(::GPU, local_rank, node_size=1)

Pin this MPI rank to one NVIDIA GPU.

`local_rank` is this process among the ranks on the same node (0, 1, 2, ...).
`node_size` is how many ranks share that node. Rank 0 uses GPU 0, rank 1 uses
GPU 1, and so on, wrapping if there are more ranks than visible GPUs.

If ranks will share a GPU, rank 0 prints a warning. The warning is skipped
when `CUDA_VISIBLE_DEVICES` already shows a single GPU.
"""
function assign_local_device!(::GPU{<:CUDABackend}, local_rank::Integer, node_size::Integer=1)
    CUDA.functional() || error(
        "MPISpec child_architecture = GPU() needs a functional CUDA device on this rank. " *
        "Load CUDA (`using CUDA`) and run on a GPU node.",
    )
    ndev = length(CUDA.devices())
    ndev < 1 && error("No CUDA devices visible to this MPI rank.")
    vis = get(ENV, "CUDA_VISIBLE_DEVICES", "")
    # `--gpus-per-task=1` will expose one device per rank, that is not sharing.
    isolated = ndev == 1 && !isempty(vis) && !occursin(',', vis)
    if node_size > ndev && local_rank == 0 && !isolated
        @warn "This node has $node_size MPI ranks and $ndev GPU(s). " *
              "Ranks will share devices (node-local rank modulo $ndev). " *
              "Request one GPU per rank, or reduce the number of ranks on the node."
    end
    CUDA.device!(Int(local_rank) % ndev)
    return nothing
end

function __init__()
    if CUDA.functional()
        WAVI.Architectures.register_gpu_backend!(GPU())
    end
end

end
