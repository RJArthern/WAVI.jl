using Parameters
using MPI

using WAVI.Parameters
using WAVI.Stencils

import WAVI: AbstractField, AbstractGrid, AbstractMeltRate, AbstractModel
import WAVI.Fields: GridField, InitialConditions, HGrid, UGrid, VGrid, CGrid, SigmaGrid
import WAVI.Grids: Grid
import WAVI.MeltRates: UniformMeltRate
import WAVI.Models: BasicSpec, Model, get_bed_elevation
import WAVI.Processes: update_state!, update_model_velocities!, update_velocities!, update_velocities_on_h_grid!
import WAVI.Wavelets: UWavelets, VWavelets

# MPI send/recv, Gatherv, and scalar halo blends use host memory on purpose
# (not CUDA-aware MPI). `_host(::Array)` is a no-copy so the CPU path is unchanged.

##
# Additional MPI functionality
#

"""
    mpi_sync_halos_after_thickness!(model)

Synchronises halo regions after a thickness update.
If Partition-of-Unity (`pou=true`) is enabled, only thickness `h` is exchanged to preserve the assembled velocity field. Otherwise, `h`, `u`, and `v` are all synchronised.
"""
function mpi_sync_halos_after_thickness!(model::AbstractModel{<:Any, <:Any, <:MPISpec})
    if model.spec.pou
        halo_exchange!(model; fields=[:h])
    else
        halo_exchange!(model)
    end
    return nothing
end

"""
    mpi_sync_halos_initial!(model)

Performs a one-time synchronization of `h`, `u`, and `v` halo regions before the first time step.
"""
function mpi_sync_halos_initial!(model::AbstractModel{<:Any, <:Any, <:MPISpec})
    halo_exchange!(model)
    return nothing
end

"""
    apply_halo_exchange_blends!(...)

Blends newly received values from neighbours into the local halo cells.
At the subgrid corners, it computes a smooth weighted mix from both adjacent neighbours.
"""
function apply_halo_exchange_blends!(
    local_field::AbstractMatrix,
    L0::AbstractMatrix,
    recv_left::Union{Nothing, AbstractMatrix},
    recv_right::Union{Nothing, AbstractMatrix},
    recv_top::Union{Nothing, AbstractMatrix},
    recv_bottom::Union{Nothing, AbstractMatrix},
    W_left::AbstractVector,
    W_right::AbstractVector,
    W_top::AbstractMatrix,
    W_bottom::AbstractMatrix,
    lh::Int,
    rh::Int,
    th::Int,
    bh::Int,
    nx::Int,
    ny::Int,
    left::Int,
    right::Int,
    top::Int,
    bottom::Int,
)
    j_mid = (th + 1):(ny - bh)
    i_mid = (lh + 1):(nx - rh)

    if left > -1 && lh > 0
        if isempty(j_mid)
            local_field[1:lh, :] .= (1.0 .- W_left) .* recv_left .+ W_left .* L0[1:lh, :]
        else
            local_field[1:lh, j_mid] .=
                (1.0 .- W_left) .* recv_left[:, j_mid] .+ W_left .* L0[1:lh, j_mid]
        end
    end
    if right > -1 && rh > 0
        ir = (nx - rh + 1):nx
        if isempty(j_mid)
            local_field[ir, :] .= (1.0 .- W_right) .* recv_right .+ W_right .* L0[ir, :]
        else
            local_field[ir, j_mid] .=
                (1.0 .- W_right) .* recv_right[:, j_mid] .+ W_right .* L0[ir, j_mid]
        end
    end
    if top > -1 && th > 0
        if isempty(i_mid)
            local_field[:, 1:th] .= (1.0 .- W_top) .* recv_top .+ W_top .* L0[:, 1:th]
        else
            local_field[i_mid, 1:th] .=
                (1.0 .- W_top) .* recv_top[i_mid, :] .+ W_top .* L0[i_mid, 1:th]
        end
    end
    if bottom > -1 && bh > 0
        jb = (ny - bh + 1):ny
        if isempty(i_mid)
            local_field[:, jb] .= (1.0 .- W_bottom) .* recv_bottom .+ W_bottom .* L0[:, jb]
        else
            local_field[i_mid, jb] .=
                (1.0 .- W_bottom) .* recv_bottom[i_mid, :] .+ W_bottom .* L0[i_mid, jb]
        end
    end

    if left > -1 && top > -1 && lh > 0 && th > 0
        for j in 1:th, i in 1:lh
            wl, wt = W_left[i], W_top[1, j]
            RL, RT = recv_left[i, j], recv_top[i, j]
            local_field[i, j] = (1 - wl) * (1 - wt) * RL + (1 - wl) * wt * RT + wl * L0[i, j]
        end
    end
    if right > -1 && top > -1 && rh > 0 && th > 0
        for j in 1:th, i_loc in 1:rh
            i = nx - rh + i_loc
            wr, wt = W_right[i_loc], W_top[1, j]
            RR, RT = recv_right[i_loc, j], recv_top[i, j]
            local_field[i, j] = (1 - wr) * (1 - wt) * RR + (1 - wr) * wt * RT + wr * L0[i, j]
        end
    end
    if left > -1 && bottom > -1 && lh > 0 && bh > 0
        for j_loc in 1:bh, i in 1:lh
            j = ny - bh + j_loc
            wl, wb = W_left[i], W_bottom[1, j_loc]
            RL, RB = recv_left[i, j], recv_bottom[i, j_loc]
            local_field[i, j] = (1 - wl) * (1 - wb) * RL + (1 - wl) * wb * RB + wl * L0[i, j]
        end
    end
    if right > -1 && bottom > -1 && rh > 0 && bh > 0
        for j_loc in 1:bh, i_loc in 1:rh
            i = nx - rh + i_loc
            j = ny - bh + j_loc
            wr, wb = W_right[i_loc], W_bottom[1, j_loc]
            RR, RB = recv_right[i_loc, j], recv_bottom[i, j_loc]
            local_field[i, j] = (1 - wr) * (1 - wb) * RR + (1 - wr) * wb * RB + wr * L0[i, j]
        end
    end
    return nothing
end

"""
    mpi_velocity_pou_weights(model::AbstractModel{<:Any,<:Any,<:MPISpec})

Return partition-of-unity (PoU) weights for the local U and V grids.

Weights are computed once and cached on `model.spec.pou_scratch` (with strip buffers
and prolong workspaces) for the life of the local velocity arrays.
"""
function mpi_velocity_pou_weights(model::AbstractModel{<:Any, <:Any, <:MPISpec})
    scratch = ensure_mpi_pou_scratch!(model)
    return scratch.ωu, scratch.ωv
end

"""
Build or reuse `MPIPoUScratch` on `model.spec` for the current local `u`/`v` sizes.
"""
function ensure_mpi_pou_scratch!(model::AbstractModel{<:Any, <:Any, <:MPISpec})
    @unpack gu, gv = model.fields
    spec = model.spec
    @unpack halo, top, right, bottom, left = spec
    mu, nu = size(gu.u)
    mv, nv = size(gv.v)
    T = eltype(gu.u)

    scratch = spec.pou_scratch
    host_u = _is_host_array(gu.u)
    if scratch !== nothing &&
       size(scratch.ωu) == (mu, nu) &&
       size(scratch.ωv) == (mv, nv) &&
       eltype(scratch.ωu) === T &&
       _is_host_array(scratch.ωu) === host_u
        return scratch
    end

    # Keep weights at 1.0 on true domain boundaries
    leavei1 = left < 0
    leaveim = right < 0
    leavej1 = top < 0
    leavejn = bottom < 0

    # Use ramp width of 2*halo to ensure weights sum to exactly 1.0 across overlapping ranks
    o = 2 * halo
    ωu = partition_of_unity(mu, nu, leavei1, leaveim, leavej1, leavejn, o, o - 1)
    ωv = partition_of_unity(mv, nv, leavei1, leaveim, leavej1, leavejn, o - 1, o)
    if host_u
        scratch = MPIPoUScratch(Matrix{T}(ωu), Matrix{T}(ωv))
    else
        ωu_dev = similar(gu.u, T, size(ωu)...)
        ωv_dev = similar(gv.v, T, size(ωv)...)
        copyto!(ωu_dev, ωu)
        copyto!(ωv_dev, ωv)
        scratch = MPIPoUScratch(ωu_dev, ωv_dev)
    end
    spec.pou_scratch = scratch
    return scratch
end

@inline function ensure_strip_buf!(buf::Vector{T}, n::Int) where {T}
    length(buf) < n && resize!(buf, n)
    return buf
end

"""
Build or reuse `MPIHaloScratch` on `model.spec` for RAS send/recv packs.
"""
function ensure_mpi_halo_scratch!(model::AbstractModel{<:Any, <:Any, <:MPISpec})
    spec = model.spec
    T = eltype(model.fields.gh.h)
    scratch = spec.halo_scratch
    if scratch !== nothing && eltype(scratch.send_l) === T && length(scratch.W_left) == spec.halo
        return scratch
    end
    scratch = MPIHaloScratch(T, spec.halo, spec.damping)
    spec.halo_scratch = scratch
    return scratch
end

"""
    copy_ras_l0!(scratch, field_data, local_field, fields)

Save a CPU copy of the current local field before neighbour values are mixed in.
Damping uses this saved copy so the update is a blend of old and new values,
not a full overwrite. Returns the saved CPU matrix. CPU fields are already on
the host, so that copy is free.
"""
function copy_ras_l0!(scratch::MPIHaloScratch{T}, field_data, local_field, fields) where {T}
    nx, ny = size(local_field)
    buf = ras_l0_workspace!(scratch, field_data, fields, nx, ny)
    copyto!(buf, _host(local_field))
    return buf
end

"""
    ras_l0_workspace!(scratch, field_data, fields, nx, ny)

Get a reusable CPU workspace for the thickness (`h`), `u`, or `v` field.
The matrix is resized only when the local patch shape changes. Both the CPU
and GPU halo paths use this workspace when mixing neighbour values.
"""
function ras_l0_workspace!(scratch::MPIHaloScratch{T}, field_data, fields, nx, ny) where {T}
    if field_data === fields.gh
        scratch.l0_h = _ensure_matrix!(scratch.l0_h, nx, ny)
        return scratch.l0_h
    elseif field_data === fields.gu
        scratch.l0_u = _ensure_matrix!(scratch.l0_u, nx, ny)
        return scratch.l0_u
    else
        scratch.l0_v = _ensure_matrix!(scratch.l0_v, nx, ny)
        return scratch.l0_v
    end
end

function _ensure_matrix!(buf::Matrix{T}, nx, ny) where {T}
    size(buf) == (nx, ny) && return buf
    return Matrix{T}(undef, nx, ny)
end

"""
    pack_halo_strip!(buf, field, irange, jrange, scratch=nothing)

Copy one edge strip of `field` into a flat CPU buffer ready for MPI.

On the CPU this is a direct copy of that strip. On the GPU the strip is packed
on the device first, then only that strip is copied to the CPU. Returns a view
of the packed values in `buf`.
"""
function pack_halo_strip!(buf::Vector{T}, field, irange, jrange, scratch=nothing) where {T}
    ni, nj = length(irange), length(jrange)
    n = ni * nj
    n == 0 && return view(buf, 1:0)
    ensure_strip_buf!(buf, n)
    if _is_host_array(field)
        copyto!(reshape(view(buf, 1:n), ni, nj), @view field[irange, jrange])
    else
        pack_halo_strip_device!(buf, field, irange, jrange, scratch)
    end
    return view(buf, 1:n)
end

"""
    ensure_dev_halo_strip!(scratch, field, n)

Get a reusable pack buffer of length `n` that lives on the same place as
`field` (CPU or GPU). Allocates a new buffer only when the old one is missing,
too short, or on the wrong place.
"""
function ensure_dev_halo_strip!(scratch, field, n::Int)
    buf = scratch.dev_halo_strip
    if buf !== nothing &&
       buf isa AbstractVector &&
       eltype(buf) === eltype(field) &&
       length(buf) >= n &&
       _is_host_array(buf) === _is_host_array(field)
        return buf
    end
    scratch.dev_halo_strip = similar(field, n)
    return scratch.dev_halo_strip
end

"""
    pack_halo_strip_device!(buf, field, irange, jrange, scratch)

Gather one edge strip from a GPU field into a flat CPU buffer.

A small kernel packs the strip on the GPU, then the packed strip is copied to
`buf`. The copy must use full arrays with start offsets; copying through a
view can fail on CUDA.
"""
function pack_halo_strip_device!(buf::Vector, field, irange, jrange, scratch)
    ni = length(irange)
    n = ni * length(jrange)
    dev = scratch === nothing ? similar(field, n) : ensure_dev_halo_strip!(scratch, field, n)
    launch!(
        _pack_halo_strip!,
        dev,
        field,
        Int(first(irange)),
        Int(first(jrange)),
        Int(ni);
        ndrange = n,
    )
    copyto!(buf, 1, dev, 1, n)
    return nothing
end

"""
    unpack_halo_strip!(field, scratch, src, irange, jrange)

Write a flat CPU strip back into one edge block of a GPU field.
If `src` is not already a packed vector of the right length, it is packed
first, then copied to the GPU and scattered into `field`.
"""
function unpack_halo_strip!(field, scratch, src, irange, jrange)
    ni, nj = length(irange), length(jrange)
    n = ni * nj
    n == 0 && return nothing
    dev = ensure_dev_halo_strip!(scratch, field, n)
    src_vec = src isa Vector && length(src) == n ? src : begin
        tmp = Vector{eltype(src)}(undef, n)
        copyto!(tmp, vec(src))
        tmp
    end
    copyto!(dev, 1, src_vec, 1, n)
    launch!(
        _unpack_halo_strip!,
        field,
        dev,
        Int(first(irange)),
        Int(first(jrange)),
        Int(ni);
        ndrange = n,
    )
    return nothing
end

"""
    add_halo_strip!(field, scratch, src, irange, jrange)

Add a flat CPU strip into one edge block of `field`.
Host fields add in place. Device fields copy the strip onto the pack
workspace and scatter-add with `_add_halo_strip!`.
"""
function add_halo_strip!(field, scratch, src, irange, jrange)
    ni, nj = length(irange), length(jrange)
    n = ni * nj
    n == 0 && return nothing
    if _is_host_array(field)
        field[irange, jrange] .+= reshape(src, ni, nj)
        return nothing
    end
    dev = ensure_dev_halo_strip!(scratch, field, n)
    src_vec = src isa Vector && length(src) == n ? src : begin
        tmp = Vector{eltype(src)}(undef, n)
        copyto!(tmp, vec(src))
        tmp
    end
    copyto!(dev, 1, src_vec, 1, n)
    launch!(
        _add_halo_strip!,
        field,
        dev,
        Int(first(irange)),
        Int(first(jrange)),
        Int(ni);
        ndrange = n,
    )
    return nothing
end

"""
    unpack_halo_rings!(field, work, scratch, lh, rh, th, bh)

Copy the left, right, top, and bottom halo edges from the CPU workspace
`work` back onto the GPU field. Edges with zero width are skipped.
"""
function unpack_halo_rings!(field, work, scratch, lh, rh, th, bh)
    nx, ny = size(work)
    lh > 0 && unpack_halo_strip!(field, scratch, view(work, 1:lh, :), 1:lh, 1:ny)
    rh > 0 && unpack_halo_strip!(field, scratch, view(work, (nx - rh + 1):nx, :), (nx - rh + 1):nx, 1:ny)
    th > 0 && unpack_halo_strip!(field, scratch, view(work, :, 1:th), 1:nx, 1:th)
    bh > 0 && unpack_halo_strip!(field, scratch, view(work, :, (ny - bh + 1):ny), 1:nx, (ny - bh + 1):ny)
    return nothing
end

"""
    copy_halo_rings_to_host!(work, field, scratch, lh, rh, th, bh)

Copy the left, right, top, and bottom halo edges from the GPU field onto the
CPU workspace `work`, so neighbour mixing can run on the CPU. Reuses the MPI
send buffers; call this only after sends have finished.
"""
function copy_halo_rings_to_host!(work, field, scratch, lh, rh, th, bh)
    nx, ny = size(work)
    if lh > 0
        n = lh * ny
        pack_halo_strip!(scratch.send_l, field, 1:lh, 1:ny, scratch)
        copyto!(view(work, 1:lh, :), reshape(view(scratch.send_l, 1:n), lh, ny))
    end
    if rh > 0
        n = rh * ny
        pack_halo_strip!(scratch.send_r, field, (nx - rh + 1):nx, 1:ny, scratch)
        copyto!(view(work, (nx - rh + 1):nx, :), reshape(view(scratch.send_r, 1:n), rh, ny))
    end
    if th > 0
        n = nx * th
        pack_halo_strip!(scratch.send_t, field, 1:nx, 1:th, scratch)
        copyto!(view(work, :, 1:th), reshape(view(scratch.send_t, 1:n), nx, th))
    end
    if bh > 0
        n = nx * bh
        pack_halo_strip!(scratch.send_b, field, 1:nx, (ny - bh + 1):ny, scratch)
        copyto!(view(work, :, (ny - bh + 1):ny), reshape(view(scratch.send_b, 1:n), nx, bh))
    end
    return nothing
end

"""
    recv_halo_strip!(buf, n)

Make sure the CPU receive buffer is large enough for `n` values, and return a
view of that length for MPI to write into.
"""
function recv_halo_strip!(buf::Vector{T}, n::Int) where {T}
    ensure_strip_buf!(buf, n)
    return view(buf, 1:n)
end

"""
    mpi_pou_add_neighbour_strips!(field, overlapi, overlapj; scratch, ...)

Additively exchange PoU contribution strips with cardinal neighbours.

`overlapi` and `overlapj` are the Partition-of-Unity (PoU) ramp widths (e.g. `2*halo` or
`2*halo+1` depending on the staggered axis). This terminology matches `partition_of_unity`.

Communication happens in the X direction first, and the values received are immediately
added to the `field`. When the Y direction is then exchanged, it propagates
any diagonal corner contributions without requiring an explicit corner message.

Send/recv packs are reused via `scratch` (`MPIPoUScratch`). `field` is the
weighted contrib on the same backend as the velocities. Overlap strips are
packed for host MPI; this is not CUDA-aware MPI.
"""
function mpi_pou_add_neighbour_strips!(
    field::AbstractMatrix{T},
    overlapi::Int,
    overlapj::Int;
    left::Int,
    right::Int,
    top::Int,
    bottom::Int,
    comm,
    tag_base::Int,
    scratch::MPIPoUScratch{T},
) where {T}
    nx, ny = size(field)
    (overlapi >= 0 && overlapj >= 0) || throw(ArgumentError("PoU strip depths must be non-negative"))
    msg = "strip depth must be strictly less than the local array size " *
          "(prefer a coarser process grid, e.g. px×1 on narrow domains)"
    (left < 0 || overlapi < nx) || throw(ArgumentError("PoU x-strip depth overlapi=$overlapi too large for nx=$nx; $msg"))
    (right < 0 || overlapi < nx) || throw(ArgumentError("PoU x-strip depth overlapi=$overlapi too large for nx=$nx; $msg"))
    (top < 0 || overlapj < ny) || throw(ArgumentError("PoU y-strip depth overlapj=$overlapj too large for ny=$ny; $msg"))
    (bottom < 0 || overlapj < ny) || throw(ArgumentError("PoU y-strip depth overlapj=$overlapj too large for ny=$ny; $msg"))

    # --- X (left / right): exchange original local contributions ---
    if overlapi > 0 && (left > -1 || right > -1)
        requests_x = MPI.RequestSet()
        n_x = overlapi * ny
        recv_left_view = recv_right_view = nothing
        if left > -1
            send_view = pack_halo_strip!(scratch.send_l, field, 1:overlapi, 1:ny, scratch)
            recv_left_view = recv_halo_strip!(scratch.recv_l, n_x)
            push!(requests_x, MPI.Isend(send_view, left, tag_base + 0, comm))
            push!(requests_x, MPI.Irecv!(recv_left_view, left, tag_base + 1, comm))
        end
        if right > -1
            send_view = pack_halo_strip!(scratch.send_r, field, (nx - overlapi + 1):nx, 1:ny, scratch)
            recv_right_view = recv_halo_strip!(scratch.recv_r, n_x)
            push!(requests_x, MPI.Isend(send_view, right, tag_base + 1, comm))
            push!(requests_x, MPI.Irecv!(recv_right_view, right, tag_base + 0, comm))
        end
        MPI.Waitall(requests_x)
        if recv_left_view !== nothing
            add_halo_strip!(field, scratch, recv_left_view, 1:overlapi, 1:ny)
        end
        if recv_right_view !== nothing
            add_halo_strip!(field, scratch, recv_right_view, (nx - overlapi + 1):nx, 1:ny)
        end
    end

    # --- Y (top / bottom): send strips after X so corners include diagonal donors ---
    if overlapj > 0 && (top > -1 || bottom > -1)
        requests_y = MPI.RequestSet()
        n_y = nx * overlapj
        recv_top_view = recv_bottom_view = nothing
        if top > -1
            send_view = pack_halo_strip!(scratch.send_t, field, 1:nx, 1:overlapj, scratch)
            recv_top_view = recv_halo_strip!(scratch.recv_t, n_y)
            push!(requests_y, MPI.Isend(send_view, top, tag_base + 2, comm))
            push!(requests_y, MPI.Irecv!(recv_top_view, top, tag_base + 3, comm))
        end
        if bottom > -1
            send_view = pack_halo_strip!(scratch.send_b, field, 1:nx, (ny - overlapj + 1):ny, scratch)
            recv_bottom_view = recv_halo_strip!(scratch.recv_b, n_y)
            push!(requests_y, MPI.Isend(send_view, bottom, tag_base + 3, comm))
            push!(requests_y, MPI.Irecv!(recv_bottom_view, bottom, tag_base + 2, comm))
        end
        MPI.Waitall(requests_y)
        if recv_top_view !== nothing
            add_halo_strip!(field, scratch, recv_top_view, 1:nx, 1:overlapj)
        end
        if recv_bottom_view !== nothing
            add_halo_strip!(field, scratch, recv_bottom_view, 1:nx, (ny - overlapj + 1):ny)
        end
    end

    return nothing
end

"""
    mpi_pou_weighted_prolong_velocities!(model, u0, v0)

Apply Additive Schwarz with Partition-of-Unity (AS-PoU) prolongation.

Each rank forms a local weighted velocity contribution, then exchanges PoU ramp
strips with cardinal neighbours so overlapping cells hold the summed blended
field. No mid-solver gather or broadcast of the global grid.
"""
function mpi_pou_weighted_prolong_velocities!(
    model::AbstractModel{<:Any, <:Any, <:MPISpec},
    u0::AbstractMatrix,
    v0::AbstractMatrix,
)
    @unpack gu, gv = model.fields
    @unpack halo, top, right, bottom, left, comm, damping = model.spec
    d = convert(eltype(gu.u), damping)
    od = one(d) - d

    scratch = ensure_mpi_pou_scratch!(model)
    ωu, ωv = scratch.ωu, scratch.ωv
    contrib_u, contrib_v = scratch.work_u, scratch.work_v

    # Weighted contrib on the same backend as the velocities. MPI still uses
    # host overlap strips (see MPIPoUScratch).
    @. contrib_u = od * ωu * gu.u
    @. contrib_v = od * ωv * gv.v
    if !iszero(d)
        @. contrib_u += d * ωu * u0
        @. contrib_v += d * ωv * v0
    end

    # Geometric patch overlap (h-overlap = 2*halo), plus one extra face on the
    # staggered axis so the strip covers the full PoU ramp for that field.
    o = 2 * halo
    mpi_pou_add_neighbour_strips!(
        contrib_u, o + 1, o;
        left, right, top, bottom, comm, tag_base = 100, scratch,
    )
    mpi_pou_add_neighbour_strips!(
        contrib_v, o, o + 1;
        left, right, top, bottom, comm, tag_base = 200, scratch,
    )

    copyto!(gu.u, contrib_u)
    copyto!(gv.v, contrib_v)
    update_velocities_on_h_grid!(model)
    return nothing
end

"""
    halo_exchange!(model; fields=[:h, :u, :v])

Exchange halo strips with cardinal MPI neighbours, then blend into the local field.

Host fields pack from views. Device fields pack the strip on the array backend
and `copyto!` only that strip through host MPI. This is not CUDA-aware MPI.
One-rank models return before packing.
"""
function halo_exchange!(model::AbstractModel{<:Any, <:Any, <:MPISpec}; fields=[:h, :u, :v])
    @unpack halo, rank = model.spec
    if halo == 0
        rank == 0 && @warn "No halo exchange to take place, returning"
        return
    end
    mpi_has_neighbours(model.spec) || return
    _halo_exchange_neighbours!(model; fields)
    return nothing
end

function _halo_exchange_neighbours!(model::AbstractModel{<:Any, <:Any, <:MPISpec}; fields=[:h, :u, :v])
    @unpack halo, comm, top, right, bottom, left, damping = model.spec
    @unpack gh, gu, gv = model.fields

    th, rh, bh, lh = get_halos(model.spec)

    # Tags for friendly neighbourhood messaging
    top_send_tag = 1
    right_send_tag = 2
    bottom_send_tag = 3
    left_send_tag = 4
    
    # Build field list based on dynamically checking which grid holds the requested fields
    exchange_pairs = Tuple{Any, Symbol}[]
    for f in fields
        hasproperty(gh, f) && push!(exchange_pairs, (gh, f))
        hasproperty(gu, f) && push!(exchange_pairs, (gu, f))
        hasproperty(gv, f) && push!(exchange_pairs, (gv, f))
    end

    # Synchronise halo regions.
    # If damping > 0, this blends the newly received neighbour values with the old local halo values.
    # If damping = 0, it simply overwrites the local halo with the neighbour's values (standard RAS).
    scratch = ensure_mpi_halo_scratch!(model)
    W_left = scratch.W_left
    W_right = scratch.W_right
    W_top = scratch.W_top
    W_bottom = scratch.W_bottom
    skip_l0 = iszero(damping)

    # Exchange requested fields
    for (field_data, attribute) in exchange_pairs
        local_field = getproperty(field_data, attribute)

        # We can only perform halo exchange on 2D arrays
        length(size(local_field)) != 2 && continue

        field_nx, field_ny = size(local_field)
        host_field = _is_host_array(local_field)

        # --- Phase 1: X-Direction Exchange (Left/Right) ---
        requests_x = MPI.RequestSet()

        # Adjust for U-staggering (nx is +1)
        # We need to skip the shared interface face to avoid 1-index shift
        off_x = (field_data === gu) ? 1 : 0

        recv_left = recv_right = nothing
        if left > -1
            ir = (lh + 1 + off_x):(lh + halo + off_x)
            n_x = halo * field_ny
            send_left_flat = pack_halo_strip!(scratch.send_l, local_field, ir, axes(local_field, 2), scratch)
            recv_left_flat = recv_halo_strip!(scratch.recv_l, n_x)
            push!(requests_x, MPI.Isend(send_left_flat, left, left_send_tag, comm))
            push!(requests_x, MPI.Irecv!(recv_left_flat, left, right_send_tag, comm))
            recv_left = reshape(recv_left_flat, halo, field_ny)
        end
        if right > -1
            ir = (field_nx - rh - halo + 1 - off_x):(field_nx - rh - off_x)
            n_x = halo * field_ny
            send_right_flat = pack_halo_strip!(scratch.send_r, local_field, ir, axes(local_field, 2), scratch)
            recv_right_flat = recv_halo_strip!(scratch.recv_r, n_x)
            push!(requests_x, MPI.Isend(send_right_flat, right, right_send_tag, comm))
            push!(requests_x, MPI.Irecv!(recv_right_flat, right, left_send_tag, comm))
            recv_right = reshape(recv_right_flat, halo, field_ny)
        end

        MPI.Waitall(requests_x)

        # --- Phase 2: Y-Direction Exchange (Top/Bottom) ---
        requests_y = MPI.RequestSet()

        off_y = (field_data === gv) ? 1 : 0

        recv_top = recv_bottom = nothing
        if top > -1
            jr = (th + 1 + off_y):(th + halo + off_y)
            n_y = field_nx * halo
            send_top_flat = pack_halo_strip!(scratch.send_t, local_field, axes(local_field, 1), jr, scratch)
            recv_top_flat = recv_halo_strip!(scratch.recv_t, n_y)
            push!(requests_y, MPI.Isend(send_top_flat, top, top_send_tag, comm))
            push!(requests_y, MPI.Irecv!(recv_top_flat, top, bottom_send_tag, comm))
            recv_top = reshape(recv_top_flat, field_nx, halo)
        end
        if bottom > -1
            jr = (field_ny - bh - halo + 1 - off_y):(field_ny - bh - off_y)
            n_y = field_nx * halo
            send_bottom_flat = pack_halo_strip!(scratch.send_b, local_field, axes(local_field, 1), jr, scratch)
            recv_bottom_flat = recv_halo_strip!(scratch.recv_b, n_y)
            push!(requests_y, MPI.Isend(send_bottom_flat, bottom, bottom_send_tag, comm))
            push!(requests_y, MPI.Irecv!(recv_bottom_flat, bottom, top_send_tag, comm))
            recv_bottom = reshape(recv_bottom_flat, field_nx, halo)
        end

        MPI.Waitall(requests_y)

        if host_field
            blend_field = local_field
            L0 = skip_l0 ? local_field : copy_ras_l0!(scratch, field_data, local_field, model.fields)
        else
            blend_field = ras_l0_workspace!(scratch, field_data, model.fields, field_nx, field_ny)
            copy_halo_rings_to_host!(blend_field, local_field, scratch, lh, rh, th, bh)
            L0 = skip_l0 ? blend_field : copy(blend_field)
        end
        apply_halo_exchange_blends!(
            blend_field,
            L0,
            recv_left,
            recv_right,
            recv_top,
            recv_bottom,
            W_left,
            W_right,
            W_top,
            W_bottom,
            lh,
            rh,
            th,
            bh,
            field_nx,
            field_ny,
            left,
            right,
            top,
            bottom,
        )
        if !host_field
            unpack_halo_rings!(local_field, blend_field, scratch, lh, rh, th, bh)
        end
    end
end

function _collect_mpi_field_one_rank!(model, path, local_field)
    global_field = model.spec.global_fields
    for path_el in path[2:end]
        global_field = getproperty(global_field, path_el)
    end
    copyto!(global_field, local_field)
    return global_field
end

function collect_mpi_field!(model::AbstractModel{T,N,S}, path::Vector{Symbol}) where {T,N,S<:MPISpec}
    @unpack comm, coords, global_size, rank = model.spec

    # Get the full field we want to collect into from the spec, and the equivalent local field on this member
    if path[1] != :global_fields
        error("$(path) should be referring to a global field, so the first symbol should be global_fields")
    end

    local_field = model.fields
    for path_el in path[2:end]
        local_field = getproperty(local_field, path_el)
    end

    # We only handle 2D fields!
    if length(size(local_field)) != 2
        error("Trying to exchange a field ",join(string.(path), ".")," that is not 2D, this is not possible")
    end
    if global_size == 1
        return _collect_mpi_field_one_rank!(model, path, local_field)
    end

    # Establish the local grid information, with full grid information available already from global_grid
    th, rh, bh, lh = get_halos(model.spec)
    x_sz, y_sz = size(local_field)
    x_start, x_end, y_start, y_end = get_bounds(model.spec)

    @debug "[$(rank+1)/$(global_size)", join(string.(path), "."), "$((x_sz, y_sz, x_start, x_end, y_start, y_end))"

    # Determine global placement for this field's core region (field-aware for staggered grids).
    grid_sym = length(path) >= 2 ? path[2] : :gh
    # Start and end indices for the global field in x
    sx = x_start + lh
    ex = x_end - rh
    # Start and end indices for the global field in y
    sy = y_start + th
    ey = y_end - bh
    # Adjust the end indices for staggered grids
    if grid_sym == :gu
        ex += 1
    elseif grid_sym == :gv
        ey += 1
    elseif grid_sym == :gc
        ex -= 1
        ey -= 1
    end

    # Send/Gather the remote copies from the other nodes into the full field.
    # We provide the local core size and positioning in the target global field.
    field_sz = MPI.Gather(((x_sz - lh - rh, y_sz - th - bh), sx, ex, sy, ey), 0, comm)

    # MPI gather needs the core cells (halo stripped) as one packed host array.
    # A 2D interior view is not packed in memory, and a GPU field is not on the host.
    sendbuf = _host(@view local_field[(1 + lh):(end - rh), (1 + th):(end - bh)])

    if rank == 0
        global_field = model.spec.global_fields # Not named correctly
        for path_el in path[2:end]
            global_field = getproperty(global_field, path_el)
        end
        # We calculate the global grid coordinates for all ranks 
        # based on the received sizes of their core domain (ie. no halo)
        count_sizes = map(x -> prod(x[1]), field_sz)
        field_type = eltype(sendbuf)
        recv_data = Vector{field_type}(undef, sum(count_sizes))
        recv_buffer = MPI.VBuffer(recv_data, count_sizes)
        @debug "[$(rank+1)/$(global_size) ", join(string.(path), "."), "] Gathering field $((1+lh, size(local_field)[1]-rh, 1+th, size(local_field)[2]-bh)) to buffer $(size(recv_data))"
        MPI.Gatherv!(sendbuf, recv_buffer, comm)

        idxer = collect(cumsum(count_sizes))

        for proc_rank in 0:(global_size-1)
            offset = proc_rank == 0 ? 0 : idxer[proc_rank]
            proc_data = recv_data[offset+1:offset + count_sizes[proc_rank+1]]
            sx, ex, sy, ey = field_sz[proc_rank+1][2:end]
            global_field[sx:ex, sy:ey] = reshape(proc_data, field_sz[proc_rank + 1][1])
        end
        MPI.Barrier(comm)
        return global_field
    else
        @debug "[$(rank+1)/$(global_size)] Sending ", join(string.(path), "."), " data"
        MPI.Gatherv!(sendbuf, nothing, comm)
        MPI.Barrier(comm)
        return nothing
    end
end
