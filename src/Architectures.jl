module Architectures

using Adapt
using LinearAlgebra: Diagonal
using LinearMaps: LinearMap
using SparseArrays: SparseMatrixCSC
using KernelAbstractions: KernelAbstractions as KA

export AbstractArchitecture, CPU, GPU
export device, synchronise, array_type, on_architecture, child_architecture, gpu_device
export architecture, zeros_on, assign_local_device!

"""
    AbstractArchitecture

Hardware label for where arrays and kernels should live.

`CPU()` is always available. `GPU()` is added by a vendor package extension
(currently CUDA, hopefully be able to extend this to AMD ROCm and Intel's oneAPI
in future). KernelAbstractions kernels stay vendor-neutral; only the array type
and backend change.
"""
abstract type AbstractArchitecture end

"""
    CPU()

Run on the host CPU. This is the default architecture.

KernelAbstractions uses Julia threads when they are available (`julia -t N`).
"""
struct CPU <: AbstractArchitecture end

"""
    GPU(device)

Run on a single GPU whose KernelAbstractions backend is `device`.

The no-argument constructor `GPU()` is defined by a vendor extension, not in
this module, so WAVI does not depend on CUDA (or any other GPU package).

After loading `CUDA`, `GPU()` builds `GPU(CUDA.CUDABackend(always_inline=true))`.
"""
struct GPU{D} <: AbstractArchitecture
    device::D
end

# Filled by vendor extensions when a functional device is present.
const LOADED_GPU_BACKENDS = AbstractArchitecture[]

device(::CPU) = KA.CPU()
device(arch::GPU) = arch.device

"""
    synchronise(arch)

Wait until work queued on `arch` has finished.
"""
synchronise(arch::AbstractArchitecture) = KA.synchronize(device(arch))

"""
    child_architecture(arch)

Architecture that owns local arrays.

For `CPU` and `GPU` this is `arch` itself. `MPISpec` returns its
`child_architecture` field (CPU or one GPU per rank).
"""
child_architecture(arch::AbstractArchitecture) = arch

array_type(::CPU) = Array

"""
    architecture(x)

Architecture that should own local arrays for `x`.

Architectures return themselves. Specs without an architecture field
(for example `BasicSpec`) return `CPU()`. `GPUSpec` returns its GPU.
`MPISpec` returns `spec.child_architecture`.
"""
architecture(arch::AbstractArchitecture) = arch
architecture(_) = CPU()

"""
    assign_local_device!(arch, local_rank, node_size=1)

Pin this process to a device for `arch`.

`CPU` is a no-op. A vendor extension overrides this for GPU so each MPI rank
on a node uses `local_rank` modulo the number of visible devices.
`local_rank` is the rank in the node-local MPI communicator
(`MPI.COMM_TYPE_SHARED`). `node_size` is that communicator's size. If it
exceeds the number of visible GPUs, a warning is issued (ranks will share
devices), unless `CUDA_VISIBLE_DEVICES` already isolates a single GPU
per rank.
"""
assign_local_device!(::AbstractArchitecture, local_rank::Integer, node_size::Integer=1) = nothing

on_architecture(::AbstractArchitecture, x::Number) = x
on_architecture(::AbstractArchitecture, ::Nothing) = nothing
on_architecture(::CPU, a::Array) = a

"""
    zeros_on(arch, T, dims...)

Allocate a zero array of type `array_type(arch){T}` with size `dims`.

`zeros(CuArray, nx, ny)` is not defined. Call this (or `similar` of an
existing device array) instead.
"""
function zeros_on(arch::AbstractArchitecture, ::Type{T}, dims::Integer...) where {T}
    AT = array_type(arch)
    return fill!(similar(AT{T}, dims...), zero(T))
end

function register_gpu_backend!(arch::AbstractArchitecture)
    push!(LOADED_GPU_BACKENDS, arch)
    return arch
end

"""
    gpu_device(; force=false)

Return first GPU architecture registered by a loaded vendor extension.
If none is registered, warn and return `CPU()`, unless `force=true`, in which
case an error is thrown.

Use this in batch scripts so the driver doesn't need to name CUDA. The vendor
package must still be loaded in the session (`using CUDA`, or an equivalent in
`startup.jl` or the job preamble). Borrowed approach from MLDataDevices of
Lux.jl.
"""
function gpu_device(; force::Bool=false)
    if !isempty(LOADED_GPU_BACKENDS)
        return first(LOADED_GPU_BACKENDS)
    end
    if force
        error(
            "No functional GPU device found. Load a vendor package first, " *
            "for example `using CUDA`.",
        )
    end
    @warn "No GPU backend loaded; falling back to CPU. " *
          "Load a vendor package (for example `using CUDA`) to enable GPU execution."
    return CPU()
end

Base.summary(::CPU) = "CPU"
Base.summary(arch::GPU) = "GPU{$(typeof(arch.device))}"

"""
    adapt_device_array(to, x)

Copy `x` onto the array type described by `to`, or leave it unchanged.

Arguments:
- `to`: Adapt target, typically an array type such as `Array` or `CuArray`.
- `x`: Value to consider. Dense `AbstractArray`s and `Diagonal`s are adapted
  so stencil diagonals can live next to the fields. Sparse matrices,
  Kronecker `LinearMap`s, and scalars stay as they are. `Ref`s are rebuilt.

TODO: Consider creating/transferring sparse operators onto the device in a
future optimisation. We really do not want a lot of host <-> device data
transfers, but, for an initial implementation, this will do.
"""
adapt_device_array(to, x::SparseMatrixCSC) = x
adapt_device_array(to, x::LinearMap) = x
adapt_device_array(to, x::Diagonal) = Diagonal(adapt_device_array(to, x.diag))
adapt_device_array(to, x::Base.RefValue) = Ref(adapt_device_array(to, x[]))
adapt_device_array(to, x::AbstractArray) = Adapt.adapt(to, x)
adapt_device_array(to, x) = x

"""
    adapt_structure_fields(to, obj)

Rebuild `obj` with `adapt_device_array(to, field)` applied to every field.

Arguments:
- `to`: Adapt target, passed through to `adapt_device_array`.
- `obj`: Struct whose fields should be adapted (grid or wavelet storage).

The returned value has the same concrete type as `obj`. Used by
`Adapt.adapt_structure` for `HGrid`, `UGrid`, `VGrid`, `CGrid`, `SigmaGrid`,
and the wavelet structs.
"""
function adapt_structure_fields(to, obj::T) where {T}
    return T(ntuple(i -> adapt_device_array(to, getfield(obj, i)), fieldcount(T))...)
end

end
