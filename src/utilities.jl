#1D Matrix operator utility functions.
spI(n) = spdiagm(n,n, 0 => ones(n))
∂1d(n,dx) = spdiagm(n,n+1,0 => -ones(n), 1 => ones(n))/dx
c(n) = spdiagm(n,n+1,0 => ones(n), 1 => ones(n))/2
χ(n) = spdiagm(n,n+2, 1 => ones(n))

"""
    opvec(model::AbstractModel,vec::AbstractVector)

Function to multiply a vector by the momentum operator.
"""
function opvec(model::AbstractModel,vec::AbstractVector)
    @unpack gh,gu,gv,gc=model.fields
    @assert length(vec)==(gu.n+gv.n)
    uspread=gu.spread*vec[1:gu.n]
    vspread=gv.spread*vec[(gu.n+1):(gu.n+gv.n)]
    extra = gh.dimplicit[]*(gu.∂x*(gu.crop*(gu.h[:].*uspread)).+gv.∂y*(gv.crop*(gv.h[:].*vspread)))
    opvecprod=
    [
     #x-component
     gu.samp*(gu.∂x'*(2gh.dneghηav[]*(2gu.∂x*uspread .+ gv.∂y*vspread)) .+
              gu.∂y'*(gc.crop*(gc.cent'*(gh.dneghηav[]*(gc.cent*(gc.crop*( gu.∂y*uspread .+ gv.∂x*vspread )))))) .+
              gu.dnegβeff[]*uspread .+ (gu.h[:].*(gu.∂x'*(extra))) )
        ;
     #y-component
     gv.samp*(gv.∂y'*(2gh.dneghηav[]*(2gv.∂y*vspread .+ gu.∂x*uspread)) .+
              gv.∂x'*(gc.crop*(gc.cent'*(gh.dneghηav[]*(gc.cent*(gc.crop*( gv.∂x*vspread .+ gu.∂y*uspread )))))) .+
              gv.dnegβeff[]*vspread .+ (gv.h[:].*(gv.∂y'*(extra))) )
    ]
    return opvecprod
end
"""
    restrictvec(model::AbstractModel,vec::AbstractVector)

Function to restrict a vector from the fine grid to the coarse grid, used in multigrid preconditioner.
"""
function restrictvec(model::AbstractModel,vec::AbstractVector)
    @unpack wu,wv,gu,gv=model.fields
    @assert length(vec)==(gu.n+gv.n)
    vecx=vec[1:gu.n]
    vecy=vec[(gu.n+1):(gu.n+gv.n)]
    restrictvec=
    [
     #x-component
     wu.samp[]*(wu.idwt'*(gu.spread*vecx))
        ;
     #y-component
     wv.samp[]*(wv.idwt'*(gv.spread*vecy))
    ]
    return restrictvec
end
"""
    prolongvec(model::AbstractModel,waveletvec::AbstractVector)

Function to prolong a vector from the coarse grid to the fine grid, used in multigrid preconditioner.
"""
function prolongvec(model::AbstractModel,waveletvec::AbstractVector)
    @unpack wu,wv,gu,gv=model.fields
    @assert length(waveletvec)==(wu.n[]+wv.n[])
    waveletvecx = waveletvec[1:wu.n[]]
    waveletvecy = waveletvec[(wu.n[]+1):(wu.n[]+wv.n[])]
    prolongvec =
    [
     #x-component
     gu.samp*(wu.idwt*(wu.spread[]*waveletvecx))
        ;
     #y-component
     gv.samp*(wv.idwt*(wv.spread[]*waveletvecy))
    ]
    return prolongvec
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

Compute the volume above floatation
"""
volume_above_floatation(h,b,params,grid) = sum(sum(height_above_floatation.(h,b,params))) .* grid.dx .* grid.dy

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



