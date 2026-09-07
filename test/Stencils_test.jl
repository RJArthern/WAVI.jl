using WAVI
using Test
using KernelAbstractions
using LinearAlgebra: mul!
using SparseArrays
using WAVI.Stencils
using WAVI.KroneckerProducts: ⊗

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

@testset "Picard stencil applies" begin
    # Unfused stencil chain. Used only to check the fused matvec.
    function chained_momentum_matvec!(y, model, x)
        gh = model.fields.gh
        gu = model.fields.gu
        gv = model.fields.gv
        gc = model.fields.gc
        s = WAVI.stencil_scratch!(model)
        T = eltype(gh.h)
        dx_inv = one(T) / model.grid.dx
        dy_inv = one(T) / model.grid.dy
        neg_dx_inv = -dx_inv
        neg_dy_inv = -dy_inv
        twoT = T(2)
        oneT = one(T)
        neg_twoT = -twoT
        neg_oneT = -oneT

        uspread = zero(gu.u)
        vspread = zero(gv.v)
        dudx = similar(gh.h)
        dvdy = similar(gh.h)
        dudy_c = similar(gh.h, T, gc.nxc, gc.nyc)
        dvdx_c = similar(dudy_c)
        r_xx = similar(gh.h)
        r_yy = similar(gh.h)
        r_xy = similar(dudy_c)
        extra = similar(gh.h)
        d_rxx_dx = similar(gu.u)
        d_rxy_dy = similar(gu.u)
        d_ryy_dy = similar(gv.v)
        d_rxy_dx = similar(gv.v)
        taubx = similar(gu.u)
        tauby = similar(gv.v)
        qx_crop = similar(gu.u)
        qy_crop = similar(gv.v)
        dqxdx = similar(gh.h)
        dqydy = similar(gh.h)
        d_extra_dx = similar(gu.u)
        d_extra_dy = similar(gv.v)
        fx = similar(gu.u)
        fy = similar(gv.v)

        launch!(_scatter!, uspread, view(x, 1:gu.ni), s.gu_inner_indices; ndrange = gu.ni, sync = false)
        launch!(_scatter!, vspread, view(x, (gu.ni + 1):(gu.ni + gv.ni)), s.gv_inner_indices; ndrange = gv.ni)
        launch!(_diff_x!, dudx, uspread, dx_inv; ndrange = size(dudx), sync = false)
        launch!(_diff_y!, dvdy, vspread, dy_inv; ndrange = size(dvdy), sync = false)
        launch!(_diff_y_staggered!, dudy_c, uspread, dy_inv; ndrange = size(dudy_c), sync = false)
        launch!(_diff_x_staggered!, dvdx_c, vspread, dx_inv; ndrange = size(dvdx_c))
        launch!(_scale_sum!, r_xx, dudx, dvdy, twoT, oneT, gh.dneghηav[].diag, neg_twoT; ndrange = length(r_xx), sync = false)
        launch!(_scale_sum!, r_yy, dudx, dvdy, oneT, twoT, gh.dneghηav[].diag, neg_twoT; ndrange = length(r_yy), sync = false)
        launch!(_masked_scale_sum!, r_xy, dudy_c, dvdx_c, gc.dneghηav[].diag, gc.mask; ndrange = length(r_xy), sync = false)
        launch!(_scale!, taubx, uspread, gu.dnegβeff[].diag, neg_oneT; ndrange = length(taubx), sync = false)
        launch!(_scale!, tauby, vspread, gv.dnegβeff[].diag, neg_oneT; ndrange = length(tauby), sync = false)
        launch!(_masked_mul!, qx_crop, gu.h, uspread, gu.mask; ndrange = size(qx_crop), sync = false)
        launch!(_masked_mul!, qy_crop, gv.h, vspread, gv.mask; ndrange = size(qy_crop))
        launch!(_diff_xT!, d_rxx_dx, r_xx, neg_dx_inv; ndrange = size(d_rxx_dx), sync = false)
        launch!(_diff_yT_staggered!, d_rxy_dy, r_xy, neg_dy_inv; ndrange = size(d_rxy_dy), sync = false)
        launch!(_diff_yT!, d_ryy_dy, r_yy, neg_dy_inv; ndrange = size(d_ryy_dy), sync = false)
        launch!(_diff_xT_staggered!, d_rxy_dx, r_xy, neg_dx_inv; ndrange = size(d_rxy_dx), sync = false)
        launch!(_diff_x!, dqxdx, qx_crop, dx_inv; ndrange = size(dqxdx), sync = false)
        launch!(_diff_y!, dqydy, qy_crop, dy_inv; ndrange = size(dqydy))
        launch!(_add_scale!, extra, dqxdx, dqydy, gh.dimplicit[].diag; ndrange = length(extra))
        launch!(_diff_xT!, d_extra_dx, extra, neg_dx_inv; ndrange = size(d_extra_dx), sync = false)
        launch!(_diff_yT!, d_extra_dy, extra, neg_dy_inv; ndrange = size(d_extra_dy))
        @. fx = d_rxx_dx + d_rxy_dy - taubx - gu.h * d_extra_dx
        @. fy = d_ryy_dy + d_rxy_dx - tauby - gv.h * d_extra_dy
        launch!(_gather!, view(y, 1:gu.ni), fx, s.gu_inner_indices; ndrange = gu.ni, sync = false)
        launch!(_gather!, view(y, (gu.ni + 1):(gu.ni + gv.ni)), fy, s.gv_inner_indices; ndrange = gv.ni)
        return y
    end

    grid = Grid(nx = 32, ny = 24, dx = 1000.0, dy = 1000.0, nσ = 4)
    model = Model(grid = grid, bed_elevation = fill(-500.0, grid.nx, grid.ny))

    WAVI.Processes.update_av_speed!(model)
    WAVI.Processes.update_shelf_strain_rate!(model)
    WAVI.Processes.update_βeff_on_uv_grids!(model)
    WAVI.Processes.update_rheological_operators!(model)
    rhs = WAVI.Processes.get_rhs(model)

    @test all(isfinite, model.fields.gh.av_speed)
    @test all(isfinite, model.fields.gh.shelf_strain_rate)
    @test length(rhs) == model.fields.gu.ni + model.fields.gv.ni
    @test all(isfinite, rhs)

    op = WAVI.Processes.get_op(model)
    x = collect(range(1.0; stop = 2.0, length = size(op, 1)))
    y = similar(x)
    mul!(y, op, x)
    @test all(isfinite, y)
    @test length(y) == length(x)

    y_chain = similar(x)
    chained_momentum_matvec!(y_chain, model, x)
    @test y ≈ y_chain rtol = 1e-8 atol = 1e-10

    WAVI.Processes.update_velocities!(model)
    @test all(isfinite, model.fields.gu.u)
    @test all(isfinite, model.fields.gv.v)

    WAVI.Wavelets.update_wavelets!(model)
    wu, wv, gu, gv = model.fields.wu, model.fields.wv, model.fields.gu, model.fields.gv
    Wxu_f = WAVI.wavelet_matrix(wu.nxuw, wu.levels, "forward")
    Wyu_f = WAVI.wavelet_matrix(wu.nyuw, wu.levels, "forward")
    Wxv_f = WAVI.wavelet_matrix(wv.nxvw, wv.levels, "forward")
    Wyv_f = WAVI.wavelet_matrix(wv.nyvw, wv.levels, "forward")
    u_dwt = reshape((Wyu_f ⊗ Wxu_f) * (gu.crop * gu.u[:]), wu.nxuw, wu.nyuw)
    v_dwt = reshape((Wyv_f ⊗ Wxv_f) * (gv.crop * gv.v[:]), wv.nxvw, wv.nyvw)
    @test wu.wavelets ≈ u_dwt rtol = 1e-12 atol = 1e-12
    @test wv.wavelets ≈ v_dwt rtol = 1e-12 atol = 1e-12
    threshold = model.solver_params.wavelet_threshold
    @test wu.mask == (abs.(u_dwt) .>= threshold)
    @test wv.mask == (abs.(v_dwt) .>= threshold)

    n_wu, n_wv = wu.n[], wv.n[]
    n_fine, n_coarse = gu.ni + gv.ni, n_wu + n_wv
    Wxu = WAVI.wavelet_matrix(wu.nxuw, wu.levels, "reverse")
    Wyu = WAVI.wavelet_matrix(wu.nyuw, wu.levels, "reverse")
    Wxv = WAVI.wavelet_matrix(wv.nxvw, wv.levels, "reverse")
    Wyv = WAVI.wavelet_matrix(wv.nyvw, wv.levels, "reverse")
    fine = collect(range(0.2; stop = 1.7, length = n_fine))
    coarse = collect(range(-0.4; stop = 0.9, length = n_coarse))

    # Pack a 2D field into the short coarse vector using index_map.
    # True cells are numbered down each column, then left to right.
    function pack_index_map(field, index_map, n)
        out = zeros(eltype(field), n)
        @inbounds for j in axes(index_map, 2), i in axes(index_map, 1)
            k = index_map[i, j]
            k > 0 && (out[k] = field[i, j])
        end
        return out
    end

    # Unpack a short coarse vector onto the 2D wavelet grid using index_map.
    function unpack_index_map(packed, index_map)
        field = zeros(eltype(packed), size(index_map))
        @inbounds for j in axes(index_map, 2), i in axes(index_map, 1)
            k = index_map[i, j]
            k > 0 && (field[i, j] = packed[k])
        end
        return field
    end

    r_ka = WAVI.get_restrict_fun(model)(zeros(n_coarse), fine)
    r_sp = zeros(n_coarse)
    u_idwtT = reshape((sparse(Wyu') ⊗ sparse(Wxu')) * (gu.spread_inner * fine[1:gu.ni]), wu.nxuw, wu.nyuw)
    v_idwtT = reshape((sparse(Wyv') ⊗ sparse(Wxv')) * (gv.spread_inner * fine[(gu.ni + 1):n_fine]), wv.nxvw, wv.nyvw)
    r_sp[1:n_wu] .= pack_index_map(u_idwtT, wu.index_map, n_wu)
    r_sp[(n_wu + 1):n_coarse] .= pack_index_map(v_idwtT, wv.index_map, n_wv)
    @test r_ka ≈ r_sp rtol = 1e-12 atol = 1e-12

    p_ka = WAVI.get_prolong_fun(model)(zeros(n_fine), coarse)
    p_sp = zeros(n_fine)
    u_spread = vec(unpack_index_map(coarse[1:n_wu], wu.index_map))
    v_spread = vec(unpack_index_map(coarse[(n_wu + 1):n_coarse], wv.index_map))
    p_sp[1:gu.ni] .= gu.samp_inner * ((Wyu ⊗ Wxu) * u_spread)
    p_sp[(gu.ni + 1):n_fine] .= gv.samp_inner * ((Wyv ⊗ Wxv) * v_spread)
    @test p_ka ≈ p_sp rtol = 1e-12 atol = 1e-12
end

@testset "Forward Haar vs Kronecker" begin
    for (nx, ny, levels) in ((8, 6, 2), (7, 5, 2), (16, 12, 3))
        Wx = WAVI.wavelet_matrix(nx, levels, "forward")
        Wy = WAVI.wavelet_matrix(ny, levels, "forward")
        dwt = Wy ⊗ Wx
        X = reshape(collect(range(0.1; stop = 1.9, length = nx * ny)), nx, ny)
        x = vec(X)
        a, b = copy(X), similar(X)
        @test vec(WAVI.Utilities.haar_dwt!(a, b, levels)) ≈ dwt * x atol = 1e-12
    end
end
