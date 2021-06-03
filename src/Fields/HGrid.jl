#Struct to hold information on h-grid, located at cell centers.
#@with_kw struct HGrid{T <: Real, N <: Integer}
#Nx::N
#Ny::N
#x0::T = 0.0
#y0::T = 0.0
#dx::T = 1.0
#dy::T = 1.0
#xx::Array{T,2}=[x0+(i-0.5)*dx for i=1:Nx, j=1:Ny]; @assert size(xx)==(Nx,Ny)
#yy::Array{T,2}=[y0+(j-0.5)*dy for i=1:Nx, j=1:Ny]; @assert size(yy)==(Nx,Ny)
#mask::Array{Bool,2} = trues(Nx,Ny); @assert size(mask)==(Nx,Ny); @assert mask == clip(mask)
#n::N = count(mask); @assert n == count(mask)
#crop::Diagonal{T,Array{T,1}} = Diagonal(float(mask[:])); @assert crop == Diagonal(float(mask[:]))
#samp::SparseMatrixCSC{T,N} = crop[mask[:],:]; @assert samp == crop[mask[:],:]
#spread::SparseMatrixCSC{T,N} = sparse(samp'); @assert spread == sparse(samp')
#b::Array{T,2} = params.bed_elevation; @assert size(b)==(Nx,Ny)
#h::Array{T,2} = zeros(Nx,Ny); @assert size(h)==(Nx,Ny)
#s::Array{T,2} = zeros(Nx,Ny); @assert size(s)==(Nx,Ny)
#dhdt::Array{T,2} = zeros(Nx,Ny); @assert size(dhdt)==(Nx,Ny)
#accumulation::Array{T,2} = zeros(Nx,Ny); @assert size(accumulation)==(Nx,Ny)
#basal_melt::Array{T,2} = zeros(Nx,Ny); @assert size(basal_melt)==(Nx,Ny)
#haf::Array{T,2} = zeros(Nx,Ny); @assert size(haf)==(Nx,Ny)
#grounded_fraction::Array{T,2} = ones(Nx,Ny); @assert size(grounded_fraction)==(Nx,Ny)
#dsdh::Array{T,2} = ones(Nx,Ny); @assert size(dsdh)==(Nx,Ny)
#shelf_strain_rate::Array{T,2} = zeros(Nx,Ny); @assert size(shelf_strain_rate)==(Nx,Ny)
#av_speed::Array{T,2} = zeros(Nx,Ny); @assert size(av_speed)==(Nx,Ny)
#u::Array{T,2} = zeros(Nx,Ny); @assert size(u)==(Nx,Ny)
#v::Array{T,2} = zeros(Nx,Ny); @assert size(v)==(Nx,Ny)
#us::Array{T,2} = zeros(Nx,Ny); @assert size(v)==(Nx,Ny)
#vs::Array{T,2} = zeros(Nx,Ny); @assert size(v)==(Nx,Ny)
#ub::Array{T,2} = zeros(Nx,Ny); @assert size(v)==(Nx,Ny)
#vb::Array{T,2} = zeros(Nx,Ny); @assert size(v)==(Nx,Ny)
#bed_speed::Array{T,2} = zeros(Nx,Ny); @assert size(bed_speed)==(Nx,Ny)
#weertman_c::Array{T,2} = zeros(Nx,Ny); @assert size(weertman_c)==(Nx,Ny)
#β::Array{T,2} = zeros(Nx,Ny); @assert size(β)==(Nx,Ny)
#βeff::Array{T,2} = zeros(Nx,Ny); @assert size(βeff)==(Nx,Ny)
#τbed::Array{T,2} = zeros(Nx,Ny); @assert size(τbed)==(Nx,Ny)
#ηav::Array{T,2}; @assert size(ηav)==(Nx,Ny)
#quad_f1::Array{T,2} = zeros(Nx,Ny); @assert size(quad_f1)==(Nx,Ny)
#quad_f2::Array{T,2} = h./(3*ηav); @assert size(quad_f2)==(Nx,Ny)
#dneghηav::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(zeros(Nx*Ny))*crop)
#dimplicit::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(zeros(Nx*Ny))*crop)
#end



struct HGrid{T <: Real, N  <: Integer}
    Nx::N
    Ny::N
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


#Hgrid constructor. Note to self, we keep Nx and Ny here to increase transparency in velocity solve
function HGrid(;
                Nx, 
                Ny,
                mask = trues(Nx,Ny),
                b,
                h = zeros(Nx,Ny),
                ηav = zeros(Nx,Ny))
    
    #construct operators
    n = count(mask);
    crop = Diagonal(float(mask[:]));
    samp = crop[mask[:],:]; 
    spread = sparse(samp');
    dneghηav = Ref(crop*Diagonal(zeros(Nx*Ny))*crop)
    dimplicit = Ref(crop*Diagonal(zeros(Nx*Ny))*crop)
     
    #construct quantities not passed
    s = zeros(Nx,Ny)
    dhdt = zeros(Nx,Ny) 
    accumulation = zeros(Nx,Ny)
    basal_melt = zeros(Nx,Ny)
    haf = zeros(Nx,Ny)
    grounded_fraction = ones(Nx,Ny)
    dsdh = ones(Nx,Ny)
    shelf_strain_rate = zeros(Nx,Ny)
    av_speed = zeros(Nx,Ny) 
    u = zeros(Nx,Ny) 
    v = zeros(Nx,Ny);
    us = zeros(Nx,Ny); 
    vs = zeros(Nx,Ny);
    ub = zeros(Nx,Ny); 
    vb= zeros(Nx,Ny);
    bed_speed = zeros(Nx,Ny)
    weertman_c = zeros(Nx,Ny)
    β = zeros(Nx,Ny)
    βeff = zeros(Nx,Ny);
    τbed = zeros(Nx,Ny);
    ηav; @assert size(ηav)==(Nx,Ny)
    quad_f1 = zeros(Nx,Ny); 
    quad_f2 = h./(3*ηav); 


    #check sizes of everything
    @assert size(mask)==(Nx,Ny); @assert mask == clip(mask)
    @assert n == count(mask)
    @assert crop == Diagonal(float(mask[:]))
    @assert samp == crop[mask[:],:]
    @assert spread == sparse(samp')
    @assert size(b)==(Nx,Ny)
    @assert size(h)==(Nx,Ny)
    @assert size(s)==(Nx,Ny)
    @assert size(dhdt)==(Nx,Ny)
    @assert size(accumulation)==(Nx,Ny)
    @assert size(basal_melt)==(Nx,Ny)
    @assert size(haf)==(Nx,Ny)
    @assert size(grounded_fraction)==(Nx,Ny)
    @assert size(dsdh)==(Nx,Ny)
    @assert size(shelf_strain_rate)==(Nx,Ny)
    @assert size(u)==(Nx,Ny)
    @assert size(v)==(Nx,Ny)
    @assert size(av_speed)==(Nx,Ny) 
    @assert size(ub)==(Nx,Ny)
    @assert size(vb)==(Nx,Ny)
    @assert size(us)==(Nx,Ny)
    @assert size(vs)==(Nx,Ny)
    @assert size(bed_speed)==(Nx,Ny)
    @assert size(weertman_c)==(Nx,Ny)
    @assert size(β)==(Nx,Ny)
    @assert size(βeff)==(Nx,Ny)
    @assert size(τbed)==(Nx,Ny)
    @assert size(quad_f1)==(Nx,Ny)
    @assert size(quad_f2)==(Nx,Ny)

return HGrid(
            Nx,
            Ny,
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