module Utilities

using InplaceOps
using LinearAlgebra
using Parameters

using WAVI: AbstractModel
using WAVI.KroneckerProducts
using KernelAbstractions: KernelAbstractions as KA
using WAVI.Stencils

export get_op_fun, get_restrict_fun, get_prolong_fun, pos_fraction, mismip_plus_bed,
    get_glx, glen_b, get_u_mask, get_v_mask, get_c_mask, clip, get_resid, get_resid!, 
    icedraft, height_above_floatation, volume_above_floatation, spI, ∂1d, c, χ,
    stencil_scratch!

#1D Matrix operator utility functions.
spI(n) = spdiagm(n,n, 0 => ones(n))
∂1d(n,dx) = spdiagm(n,n+1,0 => -ones(n), 1 => ones(n))/dx
c(n) = spdiagm(n,n+1,0 => ones(n), 1 => ones(n))/2
χ(n) = spdiagm(n,n+2, 1 => ones(n))


"""
    stencil_scratch!(model)

Return persistent scratch for the momentum operator and Picard stencil applies.
Allocated on first use and kept for the life of the model's `GridField`.
"""
function stencil_scratch!(model::AbstractModel)
    ref = model.fields.stencil_scratch
    s = ref[]
    if s !== nothing
        return s
    end
    ref[] = allocate_stencil_scratch(model)
    return ref[]
end

function allocate_stencil_scratch(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    grid = model.grid

    # Rheological diagonals, updated in-place each Picard iterate
    gh_dneghηav_diag = gh.dneghηav[].diag
    gc_dneghηav_diag = gc.dneghηav[].diag
    gu_dnegβeff_diag = gu.dnegβeff[].diag
    gv_dnegβeff_diag = gv.dnegβeff[].diag
    gh_dimplicit_diag = gh.dimplicit[].diag

    # Inner active indices for scatter/gather (replaces samp_inner / spread_inner)
    gu_inner_indices = findall(vec(gu.mask_inner))
    gv_inner_indices = findall(vec(gv.mask_inner))

    # Preallocate intermediate variables used by op_fun and Picard applies
    usampi = similar(gu.u, gu.ni)
    vsampi = similar(gv.v, gv.ni)
    uspread = similar(gu.u)
    vspread = similar(gv.v)

    dudx = similar(gh.h)
    dvdy = similar(gh.h)
    r_xx = similar(gh.h)
    r_yy = similar(gh.h)

    dudy_c = similar(gh.h, T, gc.nxc, gc.nyc)
    dvdx_c = similar(dudy_c)
    r_xy_crop_c = similar(dudy_c)

    d_rxx_dx = similar(gu.u)
    d_rxy_dy = similar(gu.u)
    d_ryy_dy = similar(gv.v)
    d_rxy_dx = similar(gv.v)

    taubx = similar(gu.u)
    tauby = similar(gv.v)

    qx_crop = similar(gu.u)
    dqxdx = similar(gh.h)
    qy_crop = similar(gv.v)
    dqydy = similar(gh.h)
    extra = similar(gh.h)
    d_extra_dx = similar(gu.u)
    d_extra_dy = similar(gv.v)
    h_d_extra_dx = similar(gu.u)
    h_d_extra_dy = similar(gv.v)

    fx = similar(gu.u)
    fy = similar(gv.v)
    fx_sampi = similar(gu.u, gu.ni)
    fy_sampi = similar(gv.v, gv.ni)

    surf_crop = similar(gh.h)
    ones_crop = similar(gh.h)
    tmpu = similar(gu.u)
    tmpv = similar(gv.v)
    tmpui = similar(gu.u, gu.ni)
    tmpvi = similar(gv.v, gv.ni)
    u_crop = similar(gu.u)
    v_crop = similar(gv.v)
    u_h = similar(gh.h)
    v_h = similar(gh.h)
    shear_c = similar(dudy_c)
    shear_h = similar(gh.h)
    β_crop = similar(gh.h)
    gf_crop = similar(gh.h)
    denu = similar(gu.u)
    denv = similar(gv.v)
    ipolgfu = zeros(T, gu.nxu, gu.nyu)
    ipolgfv = zeros(T, gv.nxv, gv.nyv)
    hη = similar(gh.h)
    hη_c = similar(gh.h, T, gc.nxc, gc.nyc)
    rhs = zeros(T, gu.ni + gv.ni)
    f1 = zeros(T, gu.ni + gv.ni)
    f2 = zeros(T, gu.ni + gv.ni)
    f3 = zeros(T, gu.ni + gv.ni)
    sui = zeros(T, gu.ni)
    hui = zeros(T, gu.ni)
    dui = zeros(T, gu.ni)
    svi = zeros(T, gv.ni)
    hvi = zeros(T, gv.ni)
    dvi = zeros(T, gv.ni)
    uvfixed = zeros(T, gu.nxu * gu.nyu + gv.nxv * gv.nyv)

    dx_inv = one(T) / grid.dx
    dy_inv = one(T) / grid.dy
    neg_dx_inv = -dx_inv
    neg_dy_inv = -dy_inv
    twoT = T(2)
    oneT = one(T)
    neg_twoT = -twoT
    neg_oneT = -oneT
    backend = KA.get_backend(gh.h)

    function op_fun!(opvecprod::AbstractVector,inputVector::AbstractVector;vecSampled::Bool=true)
        if vecSampled
            @assert length(inputVector)==(gu.ni+gv.ni)

            # Split vector into u- and v- components
            usampi .= @view inputVector[1:gu.ni]
            vsampi .= @view inputVector[(gu.ni+1):(gu.ni+gv.ni)]

            # Spread to vectors that include all grid points within rectangular domain.
            fill!(uspread, zero(T))
            fill!(vspread, zero(T))
            launch!(_scatter!, uspread, usampi, gu_inner_indices; ndrange = length(usampi), sync = false)
            launch!(_scatter!, vspread, vsampi, gv_inner_indices; ndrange = length(vsampi), sync = false)
            KA.synchronize(backend)
        else
            # Vector already includes all grid points within rectangular domain.
            @assert length(inputVector)==(gu.nxu*gu.nyu+gv.nxv*gv.nyv)
            uspread .= reshape(@view(inputVector[1:gu.nxu*gu.nyu]), gu.nxu, gu.nyu)
            vspread .= reshape(@view(inputVector[(gu.nxu*gu.nyu+1):end]), gv.nxv, gv.nyv)
        end

        # Extensional resistive stresses
        launch!(_diff_x!, dudx, uspread, dx_inv; ndrange = size(dudx), sync = false)
        launch!(_diff_y!, dvdy, vspread, dy_inv; ndrange = size(dvdy), sync = false)
        launch!(_diff_y_staggered!, dudy_c, uspread, dy_inv; ndrange = size(dudy_c), sync = false)
        launch!(_diff_x_staggered!, dvdx_c, vspread, dx_inv; ndrange = size(dvdx_c), sync = false)
        KA.synchronize(backend)

        # r_xx = -2 D (2 ∂x u + ∂y v), r_yy = -2 D (∂x u + 2 ∂y v)
        launch!(_scale_sum!, r_xx, dudx, dvdy, twoT, oneT, gh_dneghηav_diag, neg_twoT; ndrange = length(r_xx), sync = false)
        launch!(_scale_sum!, r_yy, dudx, dvdy, oneT, twoT, gh_dneghηav_diag, neg_twoT; ndrange = length(r_yy), sync = false)

        # Shearing resistive stresses (crop D (dudy + dvdx), then negate)
        launch!(_masked_scale_sum!, r_xy_crop_c, dudy_c, dvdx_c, gc_dneghηav_diag, gc.mask; ndrange = length(r_xy_crop_c), sync = false)

        # Basal drag
        launch!(_scale!, taubx, uspread, gu_dnegβeff_diag, neg_oneT; ndrange = length(taubx), sync = false)
        launch!(_scale!, tauby, vspread, gv_dnegβeff_diag, neg_oneT; ndrange = length(tauby), sync = false)

        # Extra terms arising from Schur complement of semi implicit system (Arthern et al. 2015).
        launch!(_masked_mul!, qx_crop, gu.h, uspread, gu.mask; ndrange = size(qx_crop), sync = false)
        launch!(_masked_mul!, qy_crop, gv.h, vspread, gv.mask; ndrange = size(qy_crop), sync = false)
        KA.synchronize(backend)

        # Gradients of resistive stresses
        launch!(_diff_xT!, d_rxx_dx, r_xx, neg_dx_inv; ndrange = size(d_rxx_dx), sync = false)
        launch!(_diff_yT_staggered!, d_rxy_dy, r_xy_crop_c, neg_dy_inv; ndrange = size(d_rxy_dy), sync = false)
        launch!(_diff_yT!, d_ryy_dy, r_yy, neg_dy_inv; ndrange = size(d_ryy_dy), sync = false)
        launch!(_diff_xT_staggered!, d_rxy_dx, r_xy_crop_c, neg_dx_inv; ndrange = size(d_rxy_dx), sync = false)
        launch!(_diff_x!, dqxdx, qx_crop, dx_inv; ndrange = size(dqxdx), sync = false)
        launch!(_diff_y!, dqydy, qy_crop, dy_inv; ndrange = size(dqydy), sync = false)
        KA.synchronize(backend)

        launch!(_add_scale!, extra, dqxdx, dqydy, gh_dimplicit_diag; ndrange = length(extra))
        launch!(_diff_xT!, d_extra_dx, extra, neg_dx_inv; ndrange = size(d_extra_dx), sync = false)
        launch!(_diff_yT!, d_extra_dy, extra, neg_dy_inv; ndrange = size(d_extra_dy), sync = false)
        KA.synchronize(backend)

        @. h_d_extra_dx = gu.h * d_extra_dx
        @. h_d_extra_dy = gv.h * d_extra_dy

        # Resistive forces resolved in x and y directions
        @. fx = d_rxx_dx + d_rxy_dy - taubx - h_d_extra_dx
        @. fy = d_ryy_dy + d_rxy_dx - tauby - h_d_extra_dy

        # Resistive forces sampled at valid grid points
        launch!(_gather!, fx_sampi, fx, gu_inner_indices; ndrange = length(fx_sampi), sync = false)
        launch!(_gather!, fy_sampi, fy, gv_inner_indices; ndrange = length(fy_sampi), sync = false)
        KA.synchronize(backend)

        opvecprod[1:gu.ni] .= fx_sampi
        opvecprod[(gu.ni+1):(gu.ni+gv.ni)] .= fy_sampi

        return opvecprod
    end

    return (
        op_fun! = op_fun!,
        gu_inner_indices = gu_inner_indices,
        gv_inner_indices = gv_inner_indices,
        surf_crop = surf_crop,
        ones_crop = ones_crop,
        tmpu = tmpu,
        tmpv = tmpv,
        tmpui = tmpui,
        tmpvi = tmpvi,
        u_crop = u_crop,
        v_crop = v_crop,
        u_h = u_h,
        v_h = v_h,
        dudx = dudx,
        dvdy = dvdy,
        dudy_c = dudy_c,
        dvdx_c = dvdx_c,
        shear_c = shear_c,
        shear_h = shear_h,
        β_crop = β_crop,
        gf_crop = gf_crop,
        denu = denu,
        denv = denv,
        ipolgfu = ipolgfu,
        ipolgfv = ipolgfv,
        hη = hη,
        hη_c = hη_c,
        rhs = rhs,
        f1 = f1,
        f2 = f2,
        f3 = f3,
        sui = sui,
        hui = hui,
        dui = dui,
        svi = svi,
        hvi = hvi,
        dvi = dvi,
        uvfixed = uvfixed,
    )
end

"""
    get_op_fun(model::AbstractModel)

Returns a function that multiplies a vector by the momentum operator.
"""
function get_op_fun(model::AbstractModel)
    return stencil_scratch!(model).op_fun!
end

"""
    get_restrict_fun(model::AbstractModel)

Returns a function that restricts a vector from the fine grid to the coarse grid, 
used in multigrid preconditioner.
"""
function get_restrict_fun(model::AbstractModel{T,N}) where {T,N}
    @unpack wu,wv,gu,gv=model.fields

    #Preallocate intermediate variables used by restrict_fun
    nxnyu = gu.nxu*gu.nyu
    nxnyv = gv.nxv*gv.nyv
    nxnywu = wu.nxuw*wu.nyuw
    nxnywv = wv.nxvw*wv.nyvw
    vecx :: Vector{T} = zeros(gu.ni)
    vecy :: Vector{T} = zeros(gv.ni)
    spreadvecx :: Vector{T} = zeros(nxnyu)
    spreadvecy :: Vector{T} = zeros(nxnyv)
    bigoutx :: Vector{T} = zeros(nxnywu)
    bigouty :: Vector{T} = zeros(nxnywv)
    outx :: Vector{T} = zeros(wu.n[])
    outy :: Vector{T} = zeros(wv.n[])
    restrictvec :: Vector{T} = zeros(wu.n[]+wv.n[])

    function restrict_fun!(restrictvec::AbstractVector,vec::AbstractVector)
        @assert length(vec)==(gu.ni+gv.ni)
        vecx .= @view vec[1:gu.ni]
        vecy .= @view vec[(gu.ni+1):(gu.ni+gv.ni)]
@!      spreadvecx = gu.spread_inner*vecx
@!      spreadvecy = gv.spread_inner*vecy
@!      bigoutx = wu.idwtᵀ*spreadvecx
@!      bigouty = wv.idwtᵀ*spreadvecy
@!      outx = wu.samp[]*bigoutx
@!      outy = wv.samp[]*bigouty

        restrictvec[1:wu.n[]] .=  outx
        restrictvec[(wu.n[]+1):(wu.n[]+wv.n[])] .= outy

        return restrictvec
    end

    # Return restrict_fun as a closure
    return restrict_fun!
end

"""
    get_prolong_fun(model::AbstractModel)

Returns a function that prolongs a vector from the coarse grid to the fine grid, 
used in multigrid preconditioner.
"""
function get_prolong_fun(model::AbstractModel{T,N}) where {T,N}
    @unpack wu,wv,gu,gv=model.fields

    #Preallocate intermediate variables used by prolong_fun
    nxnyu = gu.nxu*gu.nyu
    nxnyv = gv.nxv*gv.nyv
    nxnywu = wu.nxuw*wu.nyuw
    nxnywv = wv.nxvw*wv.nyvw
    waveletvecx :: Vector{T} = zeros(wu.n[])
    waveletvecy :: Vector{T} = zeros(wv.n[])
    spreadwaveletvecx :: Vector{T} = zeros(nxnywu)
    spreadwaveletvecy :: Vector{T} = zeros(nxnywv)
    bigoutx :: Vector{T} = zeros(nxnyu)
    bigouty :: Vector{T} = zeros(nxnyv)
    outx :: Vector{T} = zeros(gu.ni)
    outy :: Vector{T} = zeros(gv.ni)
    prolongvec :: Vector{T} = zeros(gu.ni+gv.ni)

    function prolong_fun!(prolongvec::AbstractVector,waveletvec::AbstractVector)

        @assert length(waveletvec)==(wu.n[]+wv.n[])

        waveletvecx .= @view waveletvec[1:wu.n[]]
        waveletvecy .= @view waveletvec[(wu.n[]+1):(wu.n[]+wv.n[])]
@!      spreadwaveletvecx = wu.spread[]*waveletvecx
@!      spreadwaveletvecy = wv.spread[]*waveletvecy
@!      bigoutx = wu.idwt*spreadwaveletvecx
@!      bigouty = wv.idwt*spreadwaveletvecy
@!      outx = gu.samp_inner*bigoutx
@!      outy = gv.samp_inner*bigouty

        prolongvec[1:gu.ni] .= outx
        prolongvec[(gu.ni+1):(gu.ni+gv.ni)] .= outy
        
        return prolongvec
    end
    
    # Return prolong_fun as a closure
    return prolong_fun!
end
"""
pos_fraction(z1;mask=mask) -> area_fraction, area_fraction_u, area_fraction_v

Return fraction of each grid cell with function z1 above zero. Uses bilinear
interpolation of values at three nearest cell centers to represent the function.
In:
   z1:             m x n array of gridded function values.
   mask:           m x n mask   1 = valid data, 0= invalid data.
Out:
   area_fraction:   m x n array showing area fraction of interpolated z1>0 on h-grid.
   area_fraction_u: m+1 x n array showing area fraction of interpolated z1>0 on u-grid.
   area_fraction_v: m x n+1 array showing area fraction of interpolated z1>0 on v-grid.
"""
function pos_fraction(z1::AbstractArray{T,2};mask=trues(size(z1))) where {T};

    m,n=size(z1);

    area_fraction=zeros(T,m,n);
    area_fraction_u=zeros(T,m+1,n);
    area_fraction_v=zeros(T,m,n+1);

    for quadrant=1:4

        area_fraction_quadrant=zeros(T,m,n);

        #Define (x_i,y_i,z_i) triples i=1,2,3 for each point
        #and its two nearest neighbours in each of two orthogonal directions.
        #Distance unit is half a grid cell so each grid cell
        #is made up of four unit-square quadrants.
        #Quadrant 1: defined by neighbours in grid East & grid North direction
        #In local coordinate system.
        #(x_1,y_1,z_1)=(0,0,z_1);
        #(x_2,y_2,z_2)=(2,0,z_2); Grid East
        #(x_3,y_3,z_3)=(0,2,z_3); Grid North
        #Triples for quadrants 2-4 are defined by neighbours in the other
        #directions.
        #Quadrant 2: Neighbours from grid North and grid West
        #Quadrant 3: Neighbours from grid West and grid South
        #Quadrant 4: Neighbours from grid South and grid East
        #If there is no valid neighbour in any direction the central value is
        #reused.
        if quadrant == 1
              z2=z1[[2:m;m],1:n];
              z2[.!mask[[2:m;m],1:n]]=z1[.!mask[[2:m;m],1:n]];
              z3=z1[1:m,[2:n;n]];
              z3[.!mask[1:m,[2:n;n]]]=z1[.!mask[1:m,[2:n;n]]];
        elseif quadrant == 2
              z2=z1[1:m,[2:n;n]];
              z2[.!mask[1:m,[2:n;n]]]=z1[.!mask[1:m,[2:n;n]]];
              z3=z1[[1;1:(m-1)],1:n];
              z3[.!mask[[1;1:(m-1)],1:n]]=z1[.!mask[[1;1:(m-1)],1:n]];
        elseif quadrant == 3
              z2=z1[[1;1:(m-1)],1:n];
              z2[.!mask[[1;1:(m-1)],1:n]]=z1[.!mask[[1;1:(m-1)],1:n]];
              z3=z1[1:m,[1;1:(n-1)]];
              z3[.!mask[1:m,[1;1:(n-1)]]]=z1[.!mask[1:m,[1;1:(n-1)]]];
        elseif quadrant == 4
              z2=z1[1:m,[1;1:(n-1)]];
              z2[.!mask[1:m,[1;1:(n-1)]]]=z1[.!mask[1:m,[1;1:(n-1)]]];
              z3=z1[[2:m;m],1:n];
              z3[.!mask[[2:m;m],1:n]]=z1[.!mask[[2:m;m],1:n]];
        else
                error("Quadrant not defined");
        end


        #Find triples where sign changes.
        sc=findall(.!((sign.(z1) .== sign.(z2)) .& (sign.(z1) .== sign.(z3))));
        nsc=length(sc);

        #Initialise array.
        areasc=zeros(T,nsc);

        #Equation (y=ax+b) of zero contour for planar elements;
        #Planar elements are defined by (x_i,y_i,z_i) triples i=1,2,3.
        #N.B. distance unit for x and y is half a grid cell.
        #Flip x and y if needed by exchanging z2 and z3, to avoid infinite a;
        flip=(abs.(z3[sc] .- z1[sc]) .< abs.(z2[sc] .- z1[sc]));
        a=zeros(T,nsc);
        b=zeros(T,nsc);
        a[.!flip].=-(z2[sc[.!flip]].-z1[sc[.!flip]])./(z3[sc[.!flip]].-z1[sc[.!flip]]);
        b[.!flip].=-2*z1[sc[.!flip]]./(z3[sc[.!flip]].-z1[sc[.!flip]]);
        a[flip].=-(z3[sc[flip]].-z1[sc[flip]])./(z2[sc[flip]].-z1[sc[flip]]);
        b[flip].=-2*z1[sc[flip]]./(z2[sc[flip]].-z1[sc[flip]]);

        #Areas of various useful triangles and trapezoids.
        a1::Vector{T}=0.5*(b.^2)./a;
        a2::Vector{T}=0.5*((1.0 .- b).^2)./a;
        a3::Vector{T}=0.5*a.+b;

        #Test which edges of unit square are intercepted by zero contour.
        test1::Vector{Int64} = (b .> 0.0) .+ (b .> 1.0);
        test2::Vector{Int64} = ((a+b) .> 0.0) .+ ((a+b) .> 1.0);

        #There are nine possible combinations ix=1:9.
        ix= 1 .+ test1 .+ 3*test2;

        #Work out area fraction of unit square for each different combination.
        areasc[ix.==1].= 0.0;
        areasc[ix.==2].=-a1[ix.==2];
        areasc[ix.==3].=a2[ix.==3].-a1[ix.==3];
        areasc[ix.==4].=a1[ix.==4].+a3[ix.==4];
        areasc[ix.==5].=a3[ix.==5];
        areasc[ix.==6].=a2[ix.==6].+a3[ix.==6];
        areasc[ix.==7].=1.0.-a2[ix.==7].+a1[ix.==7];
        areasc[ix.==8].=1.0.-a2[ix.==8];
        areasc[ix.==9].=1.0;

        #Checks to see whether this is the area above zero or below zero.
        areasc[b.<0.0].=1.0.-areasc[b.<0.0];
        areasc[z1[sc].<0.0].=1.0.-areasc[z1[sc].<0.0];

        #Compile areas for sign changes with trivial cases.
        area_fraction_quadrant[sc].=areasc;
        area_fraction_quadrant[(z1.>0.0) .& (z2.>0.0) .& (z3.>0.0)].=1.0;
        area_fraction_quadrant[(z1.<=0.0) .& (z2.<=0.0) .& (z3.<=0.0)].=0.0;

        #Don't use any quadrant from an invalid cell of the h-grid.
        area_fraction_quadrant[.!mask].=0.0;

        #Add areas for each quadrant together.
        #N.B. distance unit for x and y is half a grid cell.
        #Summing over quadrants gives the fraction of each grid cell above zero.
        area_fraction.=area_fraction.+0.25*area_fraction_quadrant;

        if quadrant == 1
            area_fraction_u.=area_fraction_u.+0.25*[zeros(T,1,n);area_fraction_quadrant];
        elseif quadrant == 2
            area_fraction_u.=area_fraction_u.+0.25*[area_fraction_quadrant;zeros(T,1,n)];
        elseif quadrant == 3
            area_fraction_u.=area_fraction_u.+0.25*[area_fraction_quadrant;zeros(T,1,n)];
        elseif quadrant == 4
            area_fraction_u.=area_fraction_u.+0.25*[zeros(T,1,n);area_fraction_quadrant];
        else
            error("Quadrant not defined");
        end


        if quadrant == 1
            area_fraction_v.=area_fraction_v.+0.25*[zeros(T,m,1) area_fraction_quadrant];
        elseif quadrant == 2
            area_fraction_v.=area_fraction_v.+0.25*[zeros(T,m,1) area_fraction_quadrant];
        elseif quadrant == 3
            area_fraction_v.=area_fraction_v.+0.25*[area_fraction_quadrant zeros(T,m,1)];
        elseif quadrant == 4
            area_fraction_v.=area_fraction_v.+0.25*[area_fraction_quadrant zeros(T,m,1)];
        else
            error("Quadrant not defined");
        end

    end

    return area_fraction, area_fraction_u, area_fraction_v

end



#MISMIP+ bed elevation
function mismip_plus_bed(x,y)
    xbar = 300000.0
    b0 = -150.0; b2 = -728.8; b4 = 343.91; b6 = -50.75
    wc = 24000.0; fc = 4000.0; dc = 500.0
    bx(x)=b0+b2*(x/xbar)^2+b4*(x/xbar)^4+b6*(x/xbar)^6
    by(y)= dc*( (1+exp(-2(y-wc)/fc))^(-1) + (1+exp(2(y+wc)/fc))^(-1) )
    b = max(bx(x) + by(y), -720.0)
    return b
end

"""
     get_glx(model)

Return the grounding line in the form x = x(y). Assumes each y-row has at least
one grid point where ice grounded and one where ice floating.
"""
function get_glx(model)
    @unpack fields, grid = model
    glmask=diff(sign.(fields.gh.haf),dims=1).==-2 #calculate where sign of height above floating passes thru zero
    glx1=grid.xxh[1:end-1,:][glmask] #x co-ordiates upstream of grounding line
    glx2=grid.xxh[2:end,:][glmask] #x co-ordinates immediately downstream
    haf1=fields.gh.haf[1:end-1,:][glmask] #Height above floatation immediately upstream of gl
    haf2=fields.gh.haf[2:end,:][glmask]
    glx=glx1+(glx2-glx1).*(zero(haf1)-haf1)./(haf2-haf1) #interpolate between grid points to find grounding line
    return glx
end




"""
    icedraft(s,h,sea_level_wrt_geoid)

Compute ice draft.
"""
icedraft(s,h,sea_level_wrt_geoid)=max(h-(s-sea_level_wrt_geoid),zero(typeof(h-(s-sea_level_wrt_geoid))))


"""
    height_above_floatation(h,b,params)

Compute height above floatation.

"""
height_above_floatation(h,b,params) = h - (params.density_ocean/params.density_ice)*(params.sea_level_wrt_geoid - b)


"""
    volume_above_floatation(h,b,params)

Compute the volume above floatation: integrated height above floatation for cells with positive height above floatation
"""
volume_above_floatation(h,b,params,grid) = sum(sum(height_above_floatation.(h,b,params)[height_above_floatation.(h,b, params) .> 0])) .* grid.dx .* grid.dy

"""
    glen_b(temperature,damage,params)

Compute stiffness parameter B in Glen flow law.

"""

function glen_b(temperature,damage,glen_a_ref, glen_n, glen_a_activation_energy, glen_temperature_ref, gas_const)
    glen_a0 = glen_a_ref*exp(+glen_a_activation_energy/(glen_temperature_ref*gas_const) )
    glen_b = (1-damage)*( glen_a0*exp(-glen_a_activation_energy/(temperature*gas_const)) )^(-1.0/glen_n)
    return glen_b
end


"""
    get_u_mask(h_mask)

Find mask of valid grid points on u-grid corresponding to a mask defined on h-grid.

"""
function get_u_mask(h_mask)
    #include all u faces next to a selected center
    (nx,ny)=size(h_mask)
    u_mask=falses(nx+1,ny)
    u_mask[1:end-1,1:end]=u_mask[1:end-1,1:end].|h_mask
    u_mask[2:end,1:end]=u_mask[2:end,1:end].|h_mask
    return u_mask
end
"""
    get_v_mask(h_mask)

Find mask of valid grid points on v-grid corresponding to a mask defined on h-grid.

"""
function get_v_mask(h_mask)
    #include all v faces next to a selected center
    (nx,ny)=size(h_mask)
    v_mask=falses(nx,ny+1)
    v_mask[1:end,1:end-1]=v_mask[1:end,1:end-1].|h_mask
    v_mask[1:end,2:end]=v_mask[1:end,2:end].|h_mask
    return v_mask
end
"""
    get_c_mask(h_mask)

Find mask of valid grid points on c-grid corresponding to a mask defined on h-grid.

"""
function get_c_mask(h_mask)
    #select cell corners with four neighbouring cell centers in h_mask
    c_mask=h_mask[1:end-1,1:end-1] .& h_mask[1:end-1,2:end] .& h_mask[2:end,1:end-1] .& h_mask[2:end,2:end]
    return c_mask
end
"""
    clip(trial_mask)

Find mask of valid grid points on h-grid corresponding to a trial mask, also defined on h-grid.
Clip any grid points from the trial mask that cannot be used in the model.

"""
function clip(trial_mask)

    good_corners=get_c_mask(trial_mask)

    #include all centers next to a selected corner
    (nx,ny)=size(trial_mask)
    mask=falses(nx,ny)
    mask[1:end-1,1:end-1]=mask[1:end-1,1:end-1].|good_corners
    mask[1:end-1,2:end]=mask[1:end-1,2:end].|good_corners
    mask[2:end,1:end-1]=mask[2:end,1:end-1].|good_corners
    mask[2:end,2:end]=mask[2:end,2:end].|good_corners

    return mask
end


"""
     get_resid(x,op,b)

Function to return residual b - op x

"""
function get_resid(x,op,b)
    resid=similar(b)
    get_resid!(resid,x,op,b)
end


"""
     get_resid!(resid,x,op,b)

In-place function to return residual b - op x

"""
function get_resid!(resid,x,op,b)
    mem_resid= @allocated   mul!(resid,op,x)
   # println("Memory allocated in get_resid is: ", mem_resid, " bytes")
    resid .= b .- resid
end

end
