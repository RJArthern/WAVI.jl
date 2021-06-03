#Struct to hold information on h-grid, located at cell centers.
#@with_kw struct HGrid{T <: Real, N <: Integer}
#nx::N
#ny::N
#x0::T = 0.0
#y0::T = 0.0
#dx::T = 1.0
#dy::T = 1.0
#xx::Array{T,2}=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]; @assert size(xx)==(nx,ny)
#yy::Array{T,2}=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]; @assert size(yy)==(nx,ny)
#mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny); @assert mask == clip(mask)
#n::N = count(mask); @assert n == count(mask)
#crop::Diagonal{T,Array{T,1}} = Diagonal(float(mask[:])); @assert crop == Diagonal(float(mask[:]))
#samp::SparseMatrixCSC{T,N} = crop[mask[:],:]; @assert samp == crop[mask[:],:]
#spread::SparseMatrixCSC{T,N} = sparse(samp'); @assert spread == sparse(samp')
#b::Array{T,2} = params.bed_elevation; @assert size(b)==(nx,ny)
#h::Array{T,2} = zeros(nx,ny); @assert size(h)==(nx,ny)
#s::Array{T,2} = zeros(nx,ny); @assert size(s)==(nx,ny)
#dhdt::Array{T,2} = zeros(nx,ny); @assert size(dhdt)==(nx,ny)
#accumulation::Array{T,2} = zeros(nx,ny); @assert size(accumulation)==(nx,ny)
#basal_melt::Array{T,2} = zeros(nx,ny); @assert size(basal_melt)==(nx,ny)
#haf::Array{T,2} = zeros(nx,ny); @assert size(haf)==(nx,ny)
#grounded_fraction::Array{T,2} = ones(nx,ny); @assert size(grounded_fraction)==(nx,ny)
#dsdh::Array{T,2} = ones(nx,ny); @assert size(dsdh)==(nx,ny)
#shelf_strain_rate::Array{T,2} = zeros(nx,ny); @assert size(shelf_strain_rate)==(nx,ny)
#av_speed::Array{T,2} = zeros(nx,ny); @assert size(av_speed)==(nx,ny)
#u::Array{T,2} = zeros(nx,ny); @assert size(u)==(nx,ny)
#v::Array{T,2} = zeros(nx,ny); @assert size(v)==(nx,ny)
#us::Array{T,2} = zeros(nx,ny); @assert size(v)==(nx,ny)
#vs::Array{T,2} = zeros(nx,ny); @assert size(v)==(nx,ny)
#ub::Array{T,2} = zeros(nx,ny); @assert size(v)==(nx,ny)
#vb::Array{T,2} = zeros(nx,ny); @assert size(v)==(nx,ny)
#bed_speed::Array{T,2} = zeros(nx,ny); @assert size(bed_speed)==(nx,ny)
#weertman_c::Array{T,2} = zeros(nx,ny); @assert size(weertman_c)==(nx,ny)
#β::Array{T,2} = zeros(nx,ny); @assert size(β)==(nx,ny)
#βeff::Array{T,2} = zeros(nx,ny); @assert size(βeff)==(nx,ny)
#τbed::Array{T,2} = zeros(nx,ny); @assert size(τbed)==(nx,ny)
#ηav::Array{T,2}; @assert size(ηav)==(nx,ny)
#quad_f1::Array{T,2} = zeros(nx,ny); @assert size(quad_f1)==(nx,ny)
#quad_f2::Array{T,2} = h./(3*ηav); @assert size(quad_f2)==(nx,ny)
#dneghηav::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(zeros(nx*ny))*crop)
#dimplicit::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(zeros(nx*ny))*crop)
#end



struct HGrid{T <: Real, N  <: Integer}
    nx::N
    ny::N
    mask::Array{Bool,2}
    n::N 
    crop::Diagonal{T,Array{T,1}}
    samp::SparseMatrixCSC{T,N} 
    spread::SparseMatrixCSC{T,N} 
    b::Array{T,2} 
    h::Array{T,2} 
    s::Array{T,2} 
    dhdt::Array{T,2} 
    accumulation::Array{T,2} 
    basal_melt::Array{T,2}
    haf::Array{T,2} 
    grounded_fraction::Array{T,2}
    dsdh::Array{T,2}
    shelf_strain_rate::Array{T,2} 
    av_speed::Array{T,2} 
    u::Array{T,2}
    v::Array{T,2} 
    us::Array{T,2}
    vs::Array{T,2}
    ub::Array{T,2}
    vb::Array{T,2}
    bed_speed::Array{T,2} 
    weertman_c::Array{T,2}
    β::Array{T,2}
    βeff::Array{T,2} 
    τbed::Array{T,2} 
    ηav::Array{T,2}
    quad_f1::Array{T,2} 
    quad_f2::Array{T,2} 
    dneghηav::Base.RefValue{Diagonal{T,Array{T,1}}}
    dimplicit::Base.RefValue{Diagonal{T,Array{T,1}}}
end


#Hgrid constructor. Note to self, we keep nx and ny here to increase transparency in velocity solve
function HGrid(;
                nx, 
                ny,
                mask = trues(nx,ny),
                b,
                h = zeros(nx,ny),
                ηav = zeros(nx,ny))
    
    #construct operators
    n = count(mask);
    crop = Diagonal(float(mask[:]));
    samp = crop[mask[:],:]; 
    spread = sparse(samp');
    dneghηav = Ref(crop*Diagonal(zeros(nx*ny))*crop)
    dimplicit = Ref(crop*Diagonal(zeros(nx*ny))*crop)
     
    #construct quantities not passed
    s = zeros(nx,ny)
    dhdt = zeros(nx,ny) 
    accumulation = zeros(nx,ny)
    basal_melt = zeros(nx,ny)
    haf = zeros(nx,ny)
    grounded_fraction = ones(nx,ny)
    dsdh = ones(nx,ny)
    shelf_strain_rate = zeros(nx,ny)
    av_speed = zeros(nx,ny) 
    u = zeros(nx,ny) 
    v = zeros(nx,ny);
    us = zeros(nx,ny); 
    vs = zeros(nx,ny);
    ub = zeros(nx,ny); 
    vb= zeros(nx,ny);
    bed_speed = zeros(nx,ny)
    weertman_c = zeros(nx,ny)
    β = zeros(nx,ny)
    βeff = zeros(nx,ny);
    τbed = zeros(nx,ny);
    ηav; @assert size(ηav)==(nx,ny)
    quad_f1 = zeros(nx,ny); 
    quad_f2 = h./(3*ηav); 


    #check sizes of everything
    @assert size(mask)==(nx,ny); @assert mask == clip(mask)
    @assert n == count(mask)
    @assert crop == Diagonal(float(mask[:]))
    @assert samp == crop[mask[:],:]
    @assert spread == sparse(samp')
    @assert size(b)==(nx,ny)
    @assert size(h)==(nx,ny)
    @assert size(s)==(nx,ny)
    @assert size(dhdt)==(nx,ny)
    @assert size(accumulation)==(nx,ny)
    @assert size(basal_melt)==(nx,ny)
    @assert size(haf)==(nx,ny)
    @assert size(grounded_fraction)==(nx,ny)
    @assert size(dsdh)==(nx,ny)
    @assert size(shelf_strain_rate)==(nx,ny)
    @assert size(u)==(nx,ny)
    @assert size(v)==(nx,ny)
    @assert size(av_speed)==(nx,ny) 
    @assert size(ub)==(nx,ny)
    @assert size(vb)==(nx,ny)
    @assert size(us)==(nx,ny)
    @assert size(vs)==(nx,ny)
    @assert size(bed_speed)==(nx,ny)
    @assert size(weertman_c)==(nx,ny)
    @assert size(β)==(nx,ny)
    @assert size(βeff)==(nx,ny)
    @assert size(τbed)==(nx,ny)
    @assert size(quad_f1)==(nx,ny)
    @assert size(quad_f2)==(nx,ny)

return HGrid(
            nx,
            ny,
            mask,
            n,
            crop,
            samp, 
            spread,
            b,
            h,
            s,
            dhdt,
            accumulation,
            basal_melt,
            haf,
            grounded_fraction,
            dsdh,
            shelf_strain_rate,
            av_speed,
            u,
            v,
            us,
            vs,
            ub,
            vb,
            bed_speed,
            weertman_c,
            β,
            βeff,
            τbed,
            ηav,
            quad_f1,
            quad_f2,
            dneghηav,
            dimplicit)
end