using Parameters
using MPI

using WAVI.Parameters

import WAVI: AbstractField, AbstractGrid, AbstractMeltRate, AbstractModel
import WAVI.Fields: GridField, InitialConditions, HGrid, UGrid, VGrid, CGrid, SigmaGrid
import WAVI.Grids: Grid
import WAVI.MeltRates: UniformMeltRate
import WAVI.Models: BasicSpec, Model, get_bed_elevation
import WAVI.Processes: update_state!, update_model_velocities!, update_velocities!, update_velocities_on_h_grid!
import WAVI.Wavelets: UWavelets, VWavelets

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
    if scratch !== nothing &&
       size(scratch.ωu) == (mu, nu) &&
       size(scratch.ωv) == (mv, nv) &&
       eltype(scratch.ωu) === T
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
    scratch = MPIPoUScratch(Matrix{T}(ωu), Matrix{T}(ωv))
    spec.pou_scratch = scratch
    return scratch
end

@inline function ensure_strip_buf!(buf::Vector{T}, n::Int) where {T}
    length(buf) < n && resize!(buf, n)
    return buf
end

"""
    mpi_pou_add_neighbour_strips!(field, overlapi, overlapj; scratch, ...)

Additively exchange PoU contribution strips with cardinal neighbours.

`overlapi` and `overlapj` are the Partition-of-Unity (PoU) ramp widths (e.g. `2*halo` or
`2*halo+1` depending on the staggered axis). This terminology matches `partition_of_unity`.

Communication happens in the X direction first, and the values received are immediately
added to the `field`. When the Y direction is then exchanged, it propagates
any diagonal corner contributions without requiring an explicit corner message.

Send/recv packs are reused via `scratch` (`MPIPoUScratch`).
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
            send_l = ensure_strip_buf!(scratch.send_l, n_x)
            recv_l = ensure_strip_buf!(scratch.recv_l, n_x)
            copyto!(send_l, 1, reshape(@view(field[1:overlapi, :]), n_x), 1, n_x)
            send_view = @view send_l[1:n_x]
            recv_left_view = @view recv_l[1:n_x]
            push!(requests_x, MPI.Isend(send_view, left, tag_base + 0, comm))
            push!(requests_x, MPI.Irecv!(recv_left_view, left, tag_base + 1, comm))
        end
        if right > -1
            send_r = ensure_strip_buf!(scratch.send_r, n_x)
            recv_r = ensure_strip_buf!(scratch.recv_r, n_x)
            copyto!(send_r, 1, reshape(@view(field[(nx - overlapi + 1):nx, :]), n_x), 1, n_x)
            send_view = @view send_r[1:n_x]
            recv_right_view = @view recv_r[1:n_x]
            push!(requests_x, MPI.Isend(send_view, right, tag_base + 1, comm))
            push!(requests_x, MPI.Irecv!(recv_right_view, right, tag_base + 0, comm))
        end
        MPI.Waitall(requests_x)
        if recv_left_view !== nothing
            field[1:overlapi, :] .+= reshape(recv_left_view, overlapi, ny)
        end
        if recv_right_view !== nothing
            field[(nx - overlapi + 1):nx, :] .+= reshape(recv_right_view, overlapi, ny)
        end
    end

    # --- Y (top / bottom): send strips after X so corners include diagonal donors ---
    if overlapj > 0 && (top > -1 || bottom > -1)
        requests_y = MPI.RequestSet()
        n_y = nx * overlapj
        recv_top_view = recv_bottom_view = nothing
        if top > -1
            send_t = ensure_strip_buf!(scratch.send_t, n_y)
            recv_t = ensure_strip_buf!(scratch.recv_t, n_y)
            copyto!(send_t, 1, reshape(@view(field[:, 1:overlapj]), n_y), 1, n_y)
            send_view = @view send_t[1:n_y]
            recv_top_view = @view recv_t[1:n_y]
            push!(requests_y, MPI.Isend(send_view, top, tag_base + 2, comm))
            push!(requests_y, MPI.Irecv!(recv_top_view, top, tag_base + 3, comm))
        end
        if bottom > -1
            send_b = ensure_strip_buf!(scratch.send_b, n_y)
            recv_b = ensure_strip_buf!(scratch.recv_b, n_y)
            copyto!(send_b, 1, reshape(@view(field[:, (ny - overlapj + 1):ny]), n_y), 1, n_y)
            send_view = @view send_b[1:n_y]
            recv_bottom_view = @view recv_b[1:n_y]
            push!(requests_y, MPI.Isend(send_view, bottom, tag_base + 3, comm))
            push!(requests_y, MPI.Irecv!(recv_bottom_view, bottom, tag_base + 2, comm))
        end
        MPI.Waitall(requests_y)
        if recv_top_view !== nothing
            field[:, 1:overlapj] .+= reshape(recv_top_view, nx, overlapj)
        end
        if recv_bottom_view !== nothing
            field[:, (ny - overlapj + 1):ny] .+= reshape(recv_bottom_view, nx, overlapj)
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
    d = oftype(gu.u[1], damping)
    od = one(d) - d

    scratch = ensure_mpi_pou_scratch!(model)
    ωu, ωv = scratch.ωu, scratch.ωv
    contrib_u, contrib_v = scratch.work_u, scratch.work_v

    # Combined local contribution (ThreadedSpec: damp*old*ω + (1-damp)*new*ω)
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

    gu.u .= contrib_u
    gv.v .= contrib_v
    update_velocities_on_h_grid!(model)
    return nothing
end

function halo_exchange!(model::AbstractModel{<:Any, <:Any, <:MPISpec}; fields::Vector{Symbol}=[:h, :u, :v])
    @unpack halo, rank, comm, top, right, bottom, left = model.spec
    @unpack gh, gu, gv = model.fields

    if halo == 0
        rank == 0 && @warn "No halo exchange to take place, returning"
        return
    end

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
    @unpack damping = model.spec
    W_left = fill(damping, halo)
    W_right = fill(damping, halo)
    W_top = fill(damping, 1, halo)
    W_bottom = fill(damping, 1, halo)

    # Exchange requested fields
    for (field_data, attribute) in exchange_pairs
        local_field = getproperty(field_data, attribute)

        # We can only perform halo exchange on 2D arrays
        length(size(local_field)) != 2 && continue

        field_nx, field_ny = size(local_field)
        L0 = copy(local_field)

        # --- Phase 1: X-Direction Exchange (Left/Right) ---
        requests_x = MPI.RequestSet()

        # Adjust for U-staggering (nx is +1)
        # We need to skip the shared interface face to avoid 1-index shift
        off_x = (field_data === gu) ? 1 : 0

        T = eltype(local_field)
        recv_left_flat = recv_right_flat = nothing
        if left > -1
            send_left = local_field[lh+1+off_x:lh+halo+off_x, :]
            send_left_flat = copy(reshape(send_left, prod(size(send_left))))
            recv_left_flat = zeros(T, prod(size(send_left)))
            push!(requests_x, MPI.Isend(send_left_flat, left, left_send_tag, comm))
            push!(requests_x, MPI.Irecv!(recv_left_flat, left, right_send_tag, comm))
        end
        if right > -1
            send_right = local_field[field_nx-rh-halo+1-off_x:field_nx-rh-off_x, :]
            send_right_flat = copy(reshape(send_right, prod(size(send_right))))
            recv_right_flat = zeros(T, prod(size(send_right)))
            push!(requests_x, MPI.Isend(send_right_flat, right, right_send_tag, comm))
            push!(requests_x, MPI.Irecv!(recv_right_flat, right, left_send_tag, comm))
        end

        MPI.Waitall(requests_x)

        recv_left = recv_left_flat === nothing ? nothing : reshape(recv_left_flat, halo, field_ny)
        recv_right = recv_right_flat === nothing ? nothing : reshape(recv_right_flat, halo, field_ny)

        # --- Phase 2: Y-Direction Exchange (Top/Bottom) ---
        requests_y = MPI.RequestSet()

        off_y = (field_data === gv) ? 1 : 0

        recv_top_flat = recv_bottom_flat = nothing
        if top > -1
            send_top = local_field[:, th+1+off_y:th+halo+off_y]
            send_top_flat = copy(reshape(send_top, prod(size(send_top))))
            recv_top_flat = zeros(T, prod(size(send_top)))
            push!(requests_y, MPI.Isend(send_top_flat, top, top_send_tag, comm))
            push!(requests_y, MPI.Irecv!(recv_top_flat, top, bottom_send_tag, comm))
        end
        if bottom > -1
            send_bottom = local_field[:, field_ny-bh-halo+1-off_y:field_ny-bh-off_y]
            send_bottom_flat = copy(reshape(send_bottom, prod(size(send_bottom))))
            recv_bottom_flat = zeros(T, prod(size(send_bottom)))
            push!(requests_y, MPI.Isend(send_bottom_flat, bottom, bottom_send_tag, comm))
            push!(requests_y, MPI.Irecv!(recv_bottom_flat, bottom, top_send_tag, comm))
        end

        MPI.Waitall(requests_y)

        recv_top = recv_top_flat === nothing ? nothing : reshape(recv_top_flat, field_nx, halo)
        recv_bottom = recv_bottom_flat === nothing ? nothing : reshape(recv_bottom_flat, field_nx, halo)

        apply_halo_exchange_blends!(
            local_field,
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
    end
end



"""Global origin (1-based) for a local 2D field on an extended patch (includes halos)."""
function mpi_global_field_origin(spec::MPISpec)::Tuple{Int, Int}
    x_start, x_end, y_start, y_end = get_bounds(spec)
    return x_start, y_start
end

"""
    mpi_add_local_patch_to_global!(model, path, local_field)

Additive gather: each rank adds its full local patch into `global_fields` at the
global indices given by `mpi_global_field_origin`. Overlapping regions sum contributions
(ThreadedSpec-style AS-PoU on the global grid).
"""
function mpi_add_local_patch_to_global!(
    model::AbstractModel{T,N,S},
    path::Vector{Symbol},
    local_field::AbstractMatrix,
) where {T,N,S<:MPISpec}
    @unpack comm, global_size, rank = model.spec
    path[1] == :global_fields || error("path must start with :global_fields")
    global_field = model.spec.global_fields
    for p in path[2:end]
        global_field = getproperty(global_field, p)
    end
    gx0, gy0 = mpi_global_field_origin(model.spec)
    nx, ny = size(local_field)
    flat = vec(copy(local_field))
    meta = MPI.Gather((nx, ny, gx0, gy0, length(flat)), 0, comm)
    if rank == 0
        counts = [m[5] for m in meta]
        recv_buf = Vector{eltype(local_field)}(undef, sum(counts))
        MPI.Gatherv!(flat, MPI.VBuffer(recv_buf, counts), comm)
        offset = 0
        for (i, (pnx, pny, pgx0, pgy0, _)) in enumerate(meta)
            n = counts[i]
            block = recv_buf[offset+1:offset+n]
            offset += n
            global_field[pgx0:(pgx0+pnx-1), pgy0:(pgy0+pny-1)] .+= reshape(block, pnx, pny)
        end
    else
        MPI.Gatherv!(flat, nothing, comm)
    end
    MPI.Barrier(comm)
    return global_field
end

"""
    mpi_fill_local_from_global!(model, path, local_field)

Copy each rank's patch from the assembled field on `global_fields` (rank 0).

The full global array is broadcast from rank 0; each rank then indexes its own
`mpi_global_field_origin` slice.
"""
function mpi_fill_local_from_global!(
    model::AbstractModel{T,N,S},
    path::Vector{Symbol},
    local_field::AbstractMatrix,
) where {T,N,S<:MPISpec}
    @unpack comm, global_fields = model.spec
    path[1] == :global_fields || error("path must start with :global_fields")
    global_field = global_fields
    for p in path[2:end]
        global_field = getproperty(global_field, p)
    end
    MPI.Bcast!(global_field, comm)
    gx0, gy0 = mpi_global_field_origin(model.spec)
    nx, ny = size(local_field)
    local_field .= global_field[gx0:(gx0+nx-1), gy0:(gy0+ny-1)]
    return local_field
end

"""
    mpi_zero_global_field!(model, path)

Zero a global field on rank 0.
"""
function mpi_zero_global_field!(model::AbstractModel{T,N,S}, path::Vector{Symbol}) where {T,N,S<:MPISpec}
    if model.spec.rank == 0
        gf = model.spec.global_fields
        for p in path[2:end]
            gf = getproperty(gf, p)
        end
        gf .= zero(eltype(gf))
    end
    MPI.Barrier(model.spec.comm)
    return nothing
end

"""Gather disjoint core cells into `global_fields` (replace), optionally scaled."""
function mpi_init_global_core_field!(
    model::AbstractModel{T,N,S},
    path::Vector{Symbol},
    local_field::AbstractMatrix;
    scale::Real = 1,
) where {T,N,S<:MPISpec}
    @unpack comm, global_fields, global_size, rank = model.spec
    path[1] == :global_fields || error("path must start with :global_fields")
    global_field = global_fields
    for p in path[2:end]
        global_field = getproperty(global_field, p)
    end
    th, rh, bh, lh = get_halos(model.spec)
    x_sz, y_sz = size(local_field)
    x_start, x_end, y_start, y_end = get_bounds(model.spec)
    grid_sym = length(path) >= 2 ? path[2] : :gh
    sx = x_start + lh
    ex = x_end - rh
    sy = y_start + th
    ey = y_end - bh
    if grid_sym == :gu
        ex += 1
    elseif grid_sym == :gv
        ey += 1
    end
    field_sz = MPI.Gather(((x_sz - lh - rh, y_sz - th - bh), sx, ex, sy, ey), 0, comm)
    if rank == 0
        count_sizes = map(x -> prod(x[1]), field_sz)
        recv_data = Vector{eltype(local_field)}(undef, sum(count_sizes))
        recv_buffer = MPI.VBuffer(recv_data, count_sizes)
        MPI.Gatherv!(local_field[1+lh:end-rh, 1+th:end-bh], recv_buffer, comm)
        idxer = collect(cumsum(count_sizes))
        for proc_rank in 0:(global_size-1)
            offset = proc_rank == 0 ? 0 : idxer[proc_rank]
            proc_data = recv_data[offset+1:offset + count_sizes[proc_rank+1]]
            sx_p, ex_p, sy_p, ey_p = field_sz[proc_rank+1][2:end]
            global_field[sx_p:ex_p, sy_p:ey_p] .= scale .* reshape(proc_data, field_sz[proc_rank + 1][1])
        end
    else
        MPI.Gatherv!(local_field[1+lh:end-rh, 1+th:end-bh], nothing, comm)
    end
    MPI.Barrier(comm)
    return nothing
end

function collect_mpi_field!(model::AbstractModel{T,N,S}, path::Vector{Symbol}) where {T,N,S<:MPISpec}
    @unpack comm, coords, global_fields, global_size, rank = model.spec

    # Get the full field we want to collect into from the spec, and the equivalent local field on this member
    if path[1] != :global_fields
        error("$(path) should be referring to a global field, so the first symbol should be global_fields")
    end

    global_field = global_fields # Not named correctly
    local_field = model.fields
    for path_el in path[2:end] 
        global_field = getproperty(global_field, path_el)
        local_field = getproperty(local_field, path_el)
    end

    # Establish the local grid information, with full grid information available already from global_grid
    th, rh, bh, lh = get_halos(model.spec)
    # We only handle 2D fields!
    if length(size(local_field)) != 2
        error("Trying to exchange a field ",join(string.(path), ".")," that is not 2D, this is not possible")
    end
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
    
    if rank == 0
        # We calculate the global grid coordinates for all ranks 
        # based on the received sizes of their core domain (ie. no halo)
        count_sizes = map(x -> prod(x[1]), field_sz)
        field_type = eltype(local_field)
        recv_data = Vector{field_type}(undef, sum(count_sizes))
        recv_buffer = MPI.VBuffer(recv_data, count_sizes)
        @debug "[$(rank+1)/$(global_size) ", join(string.(path), "."), "] Gathering field $((1+lh, size(local_field)[1]-rh, 1+th, size(local_field)[2]-bh)) to buffer $(size(recv_data))"
        MPI.Gatherv!(local_field[1+lh:end-rh, 1+th:end-bh], recv_buffer, comm)

        idxer = collect(cumsum(count_sizes))

        for proc_rank in 0:(global_size-1)
            offset = proc_rank == 0 ? 0 : idxer[proc_rank]
            proc_data = recv_data[offset+1:offset + count_sizes[proc_rank+1]]
            sx, ex, sy, ey = field_sz[proc_rank+1][2:end]
            global_field[sx:ex, sy:ey] = reshape(proc_data, field_sz[proc_rank + 1][1])
        end
    else
        @debug "[$(rank+1)/$(global_size)] Sending ", join(string.(path), "."), " data"
        MPI.Gatherv!(local_field[1+lh:end-rh, 1+th:end-bh], nothing, comm)
    end

    MPI.Barrier(comm)
    return global_field
end
