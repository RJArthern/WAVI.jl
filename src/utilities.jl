#1D Matrix operator utility functions.
spI(n) = spdiagm(n,n, 0 => ones(n))
∂1d(n,dx) = spdiagm(n,n+1,0 => -ones(n), 1 => ones(n))/dx
c(n) = spdiagm(n,n+1,0 => ones(n), 1 => ones(n))/2
χ(n) = spdiagm(n,n+2, 1 => ones(n))

"""
    get_op_fun(model::AbstractModel)

Returns a function that multiplies a vector by the momentum operator.
"""
function get_op_fun(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    
    #Preallocate intermediate variables used by op_fun
    nxnyh :: N = gh.nxh*gh.nyh
    nxnyu :: N = gu.nxu*gu.nyu
    nxnyv :: N = gv.nxv*gv.nyv
    nxnyc :: N = gc.nxc*gc.nyc
    usamp :: Vector{T} = zeros(gu.n);                               @assert length(usamp) == gu.n
    vsamp :: Vector{T} = zeros(gv.n);                               @assert length(vsamp) == gv.n
    uspread :: Vector{T} = zeros(nxnyu);                            @assert length(uspread) == nxnyu
    vspread :: Vector{T} = zeros(nxnyv);                            @assert length(vspread) == nxnyv
    dudx :: Vector{T} = zeros(nxnyh);                               @assert length(dudx) == nxnyh
    dvdy :: Vector{T} = zeros(nxnyh);                               @assert length(dvdy) == nxnyh
    r_xx_strain_rate_sum :: Vector{T} = zeros(nxnyh);               @assert length(r_xx_strain_rate_sum) == nxnyh
    r_yy_strain_rate_sum :: Vector{T} = zeros(nxnyh);               @assert length(r_yy_strain_rate_sum) == nxnyh
    r_xx :: Vector{T} = zeros(nxnyh);                               @assert length(r_xx) == nxnyh
    r_yy :: Vector{T} = zeros(nxnyh);                               @assert length(r_yy) == nxnyh
    dudy_c :: Vector{T} = zeros(nxnyc);                             @assert length(dudy_c) == nxnyc
    dvdx_c :: Vector{T} = zeros(nxnyc);                             @assert length(dvdx_c) == nxnyc
    r_xy_strain_rate_sum_c :: Vector{T} = zeros(nxnyc);             @assert length(r_xy_strain_rate_sum_c) == nxnyc
    r_xy_strain_rate_sum_crop_c :: Vector{T} = zeros(nxnyc);        @assert length(r_xy_strain_rate_sum_crop_c) == nxnyc
    r_xy_strain_rate_sum :: Vector{T} = zeros(nxnyh);               @assert length(r_xy_strain_rate_sum) == nxnyh
    r_xy :: Vector{T} = zeros(nxnyh);                               @assert length(r_xy) == nxnyh
    r_xy_c :: Vector{T} = zeros(nxnyc);                             @assert length(r_xy_c) == nxnyc
    r_xy_crop_c :: Vector{T} = zeros(nxnyc);                        @assert length(r_xy_crop_c) == nxnyc
    d_rxx_dx :: Vector{T} = zeros(nxnyu);                           @assert length(d_rxx_dx) == nxnyu
    d_rxy_dy :: Vector{T} = zeros(nxnyu);                           @assert length(d_rxy_dy) == nxnyu
    d_ryy_dy :: Vector{T} = zeros(nxnyv);                           @assert length(d_ryy_dy) == nxnyv
    d_rxy_dx :: Vector{T} = zeros(nxnyv);                           @assert length(d_rxy_dx) == nxnyv
    taubx :: Vector{T} = zeros(nxnyu);                              @assert length(taubx) == nxnyu
    tauby :: Vector{T} = zeros(nxnyv);                              @assert length(tauby) == nxnyv
    qx :: Vector{T} = zeros(nxnyu);                                 @assert length(qx) == nxnyu
    qx_crop :: Vector{T} = zeros(nxnyu);                            @assert length(qx_crop) == nxnyu
    dqxdx :: Vector{T} = zeros(nxnyh);                              @assert length(dqxdx) == nxnyh
    qy :: Vector{T} = zeros(nxnyv);                                 @assert length(qy) == nxnyv
    qy_crop :: Vector{T} = zeros(nxnyv);                            @assert length(qy_crop) == nxnyv
    dqydy :: Vector{T} = zeros(nxnyh);                              @assert length(dqydy) == nxnyh
    divq :: Vector{T} = zeros(nxnyh);                               @assert length(divq) == nxnyh
    extra :: Vector{T} = zeros(nxnyh);                              @assert length(extra) == nxnyh
    d_extra_dx :: Vector{T} = zeros(nxnyu);                         @assert length(d_extra_dx) == nxnyu
    d_extra_dy :: Vector{T} = zeros(nxnyv);                         @assert length(d_extra_dy) == nxnyv
    h_d_extra_dx :: Vector{T} = zeros(nxnyu);                       @assert length(h_d_extra_dx) == nxnyu
    h_d_extra_dy :: Vector{T} = zeros(nxnyv);                       @assert length(h_d_extra_dy) == nxnyv
    fx :: Vector{T} = zeros(nxnyu);                                 @assert length(fx) == nxnyu
    fy :: Vector{T} = zeros(nxnyv);                                 @assert length(fy) == nxnyv
    fx_samp :: Vector{T} = zeros(gu.n);                             @assert length(fx_samp) == gu.n
    fy_samp :: Vector{T} = zeros(gv.n);                             @assert length(fy_samp) == gv.n
    opvecprod :: Vector{T} = zeros(gu.n+gv.n);                      @assert length(opvecprod) == gu.n + gv.n

    function op_fun(vec::AbstractVector)

            @assert length(vec)==(gu.n+gv.n)

            #Split vector into u- and v- components
            usamp .= @view vec[1:gu.n]
            vsamp .= @view vec[(gu.n+1):(gu.n+gv.n)]

            #Spread to vectors that include all grid points within rectangular domain.
        @!  uspread = gu.spread*usamp
        @!  vspread = gv.spread*vsamp

            #Extensional resistive stresses
        @!  dudx = gu.∂x*uspread
        @!  dvdy = gv.∂y*vspread
            r_xx_strain_rate_sum .= 2dudx .+ dvdy
            r_yy_strain_rate_sum .= 2dvdy .+ dudx
        @!  r_xx = -2gh.dneghηav[]*r_xx_strain_rate_sum
        @!  r_yy = -2gh.dneghηav[]*r_yy_strain_rate_sum
            
            #Shearing resistive stresses
        @!  dudy_c = gu.∂y*uspread
        @!  dvdx_c = gv.∂x*vspread
            r_xy_strain_rate_sum_c .= dudy_c .+ dvdx_c
        @!  r_xy_strain_rate_sum_crop_c = gc.crop*r_xy_strain_rate_sum_c
        @!  r_xy_strain_rate_sum = gc.cent*r_xy_strain_rate_sum_crop_c
        @!  r_xy = -gh.dneghηav[]*r_xy_strain_rate_sum
        @!  r_xy_c = gc.cent'*r_xy
        @!  r_xy_crop_c = gc.crop*r_xy_c

            #Gradients of resisitve stresses
        @!  d_rxx_dx = -gu.∂x'*r_xx
        @!  d_rxy_dy = -gu.∂y'*r_xy_crop_c
        @!  d_ryy_dy = -gv.∂y'*r_yy
        @!  d_rxy_dx = -gv.∂x'*r_xy_crop_c

            #Basal drag
        @!  taubx = -gu.dnegβeff[]*uspread
        @!  tauby = -gv.dnegβeff[]*vspread
            
            #Extra terms arising from Schur complement of semi implicit system (Arthern et al. 2015).
            qx .= gu.h[:].*uspread
        @!  qx_crop = gu.crop*qx
        @!  dqxdx = gu.∂x*qx_crop
            qy .= gv.h[:].*vspread
        @!  qy_crop = gv.crop*qy
        @!  dqydy = gv.∂y*qy_crop
            divq .= dqxdx .+ dqydy
        @!  extra = gh.dimplicit[]*divq
        @!  d_extra_dx = -gu.∂x'*extra
        @!  d_extra_dy = -gv.∂y'*extra
            h_d_extra_dx .= gu.h[:].*d_extra_dx
            h_d_extra_dy .= gv.h[:].*d_extra_dy

            #Resistive forces resolved in x anf y directions
            fx .= d_rxx_dx .+ d_rxy_dy .- taubx .- h_d_extra_dx
            fy .= d_ryy_dy .+ d_rxy_dx .- tauby .- h_d_extra_dy

            #Resistive forces sampled at valid grid points
        @!  fx_samp = gu.samp*fx
        @!  fy_samp = gv.samp*fy

            opvecprod .=
            [
            #x-component
            fx_samp
                ;
            #y-component
            fy_samp
            ]
            return opvecprod
    end

    #Return op_fun as a closure
    return op_fun
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
    vecx :: Vector{T} = zeros(gu.n)
    vecy :: Vector{T} = zeros(gv.n)
    spreadvecx :: Vector{T} = zeros(nxnyu)
    spreadvecy :: Vector{T} = zeros(nxnyv)
    bigoutx :: Vector{T} = zeros(nxnywu)
    bigouty :: Vector{T} = zeros(nxnywv)
    outx :: Vector{T} = zeros(wu.n[])
    outy :: Vector{T} = zeros(wv.n[])
    restrictvec :: Vector{T} = zeros(wu.n[]+wv.n[])

    function restrict_fun(vec::AbstractVector)
        @assert length(vec)==(gu.n+gv.n)
        vecx .= @view vec[1:gu.n]
        vecy .= @view vec[(gu.n+1):(gu.n+gv.n)]
@!      spreadvecx = gu.spread*vecx
@!      spreadvecy = gv.spread*vecy
@!      bigoutx = wu.idwt'*spreadvecx
@!      bigouty = wv.idwt'*spreadvecy
@!      outx = wu.samp[]*bigoutx
@!      outy = wv.samp[]*bigouty

     
        restrictvec .=
        [
        #x-component
        outx
            ;
        #y-component
        outy
        ]
        return restrictvec
    end

    # Return restrict_fun as a closure
    return restrict_fun
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
    outx :: Vector{T} = zeros(gu.n)
    outy :: Vector{T} = zeros(gv.n)
    prolongvec :: Vector{T} = zeros(gu.n+gv.n)

    function prolong_fun(waveletvec::AbstractVector)

        @assert length(waveletvec)==(wu.n[]+wv.n[])

        waveletvecx .= @view waveletvec[1:wu.n[]]
        waveletvecy .= @view waveletvec[(wu.n[]+1):(wu.n[]+wv.n[])]
@!      spreadwaveletvecx = wu.spread[]*waveletvecx
@!      spreadwaveletvecy = wv.spread[]*waveletvecy
@!      bigoutx = wu.idwt*spreadwaveletvecx
@!      bigouty = wv.idwt*spreadwaveletvecy
@!      outx = gu.samp*bigoutx
@!      outy = gv.samp*bigouty

        prolongvec .=
        [
        #x-component
        outx
            ;
        #y-component
        outy
        ]
        return prolongvec
    end
    
    # Return prolong_fun as a closure
    return prolong_fun
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



