module Architectures

using Adapt
using LinearAlgebra: Diagonal
using LinearMaps: LinearMap
using SparseArrays: SparseMatrixCSC
using KernelAbstractions: KernelAbstractions as KA

export AbstractArchitecture, CPU, GPU
export device, synchronise, array_type, on_architecture, child_architecture, gpu_device

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

Architecture that owns local arrays. For `CPU` and `GPU` this is `arch` itself.
Aim is for distributed specs to override this in the future (e.g. one GPU per
MPI rank).
"""
child_architecture(arch::AbstractArchitecture) = arch

array_type(::CPU) = Array

on_architecture(::AbstractArchitecture, x::Number) = x
on_architecture(::AbstractArchitecture, ::Nothing) = nothing
on_architecture(::CPU, a::Array) = a

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
- `x`: Value to consider. Dense `AbstractArray`s are adapted. Sparse matrices,
  `Diagonal`s, Kronecker `LinearMap`s, `Ref`s, and scalars stay as they are so
  host operators are not copied onto the device.

TODO: Consider copying operators onto the device in a future optimisation.
We really do not want a lot of host <-> device data transfers, but, for an
initial implementation, this will do.
"""
adapt_device_array(to, x::SparseMatrixCSC) = x
adapt_device_array(to, x::Diagonal) = x
adapt_device_array(to, x::LinearMap) = x
adapt_device_array(to, x::Base.RefValue) = x
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
