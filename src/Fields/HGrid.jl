struct HGrid{T <: Real, N  <: Integer}
                   Nx :: N                                     # Number of grid cells in x-direction in HGrid
                   Ny :: N                                     # Number of grid cells in y-direction in HGrid
                 mask :: Array{Bool,2}                         # Mask specifying the model domain
                    n :: N                                     # Total number of cells in the model domain
                 crop :: Diagonal{T,Array{T,1}}                # Crop matrix: diagonal matrix with mask entries on diag
                 samp :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to model domain 
               spread :: SparseMatrixCSC{T,N}                  # Sparse form of the sampling matrix 
                    b :: Array{T,2}                            # Bed elevation
                    h :: Array{T,2}                            # Ice thickness 
                    s :: Array{T,2}                            # Current surface elevation
                 dhdt :: Array{T,2}                            # Time rate of change of thickness 
         accumulation :: Array{T,2}                            # Accumulation rate 
           basal_melt :: Array{T,2}                            # Basal melt rate    
                  haf :: Array{T,2}                            # Grid cell height above floatation
    grounded_fraction :: Array{T,2}                            # Grid cell grounded fraction 
                 dsdh :: Array{T,2}                            # Cange of surface elevation per unit thickness change
    shelf_strain_rate :: Array{T,2}                            # Strain rate appropriate for shelf (no basal drag) 
             av_speed :: Array{T,2}                            # Depth averaged speed 
                    u :: Array{T,2}                            # Depth averaged x-velocity 
                    v :: Array{T,2}                            # Depth averaged y-velocity 
                   us :: Array{T,2}                            # x-velocity at the surface 
                   vs :: Array{T,2}                            # y-velocity at the surface 
                   ub :: Array{T,2}                            # x-velocity at the bed 
                   vb :: Array{T,2}                            # y-velocity at the bed
            bed_speed :: Array{T,2}                            # Ice speed at the bed
           weertman_c :: Array{T,2}                            # Weertman drag coefficients 
                    β :: Array{T,2}                            # Raw β value (eqn 8 in Arthern 2015 JGeophysRes)
                 βeff :: Array{T,2}                            # Effective β value (eqn 12 in Arthern 2015 JGeophysRes)
                 τbed :: Array{T,2}                            # Stress at the bed
                  ηav :: Array{T,2}                            # Depth averaged viscosity
              quad_f1 :: Array{T,2}                            # F1 quadratrue field (eqn 7 in Arthern 2015 JGeophysRes)
              quad_f2 :: Array{T,2}                            # F2 quadrature field (eqn 7 in Arthern 2015 JGeophysRes)
             dneghηav :: Base.RefValue{Diagonal{T,Array{T,1}}} # Rheological operator (-h × ηav)
            dimplicit :: Base.RefValue{Diagonal{T,Array{T,1}}} # Rheological operator (-ρi × g × dt × dshs)
end


"""
    HGrid(;
            Nx, 
            Ny,
            mask = trues(Nx,Ny),
            b,
            h,
            ηav = zeros(Nx,Ny))

Construct a WAVI.jl HGrid with size (Nx,Ny)
HGrid stores fields that are defined on the problem's H grid. 
(Co-ordinates of HGrid stored in a Grid under xxh, yyh fields)

Keyword arguments
=================
    - 'Nx': (required) Number of grid cells in x-direction in HGrid (should be same as grid.nx)
            Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
    - 'Ny': (required) Number of grid cells in y-direction in HGrid (should be same as grid.ny)
    - 'mask': Mask specifying the model domain
    - 'b': (requried) Bed elevation (bottom bathymetry)
    - 'h': (required) initial thickness of the ice
    - 'ηav': depth averaged visosity initially
"""


function HGrid(;
                Nx, 
                Ny,
                mask = trues(Nx,Ny),
                b,
                h = zeros(Nx,Ny),
                ηav = zeros(Nx,Ny))

    #check the sizes of inputs
    (size(mask) == size(b) == size(h) == size(ηav) == (Nx,Ny)) || throw(DimensionMismatch("Sizes of inputs to UGrid must all be equal to Nx x Ny (i.e. $Nx x $Ny)"))

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
    v = zeros(Nx,Ny)
    us = zeros(Nx,Ny) 
    vs = zeros(Nx,Ny)
    ub = zeros(Nx,Ny) 
    vb= zeros(Nx,Ny)
    bed_speed = zeros(Nx,Ny)
    weertman_c = zeros(Nx,Ny)
    β = zeros(Nx,Ny)
    βeff = zeros(Nx,Ny)
    τbed = zeros(Nx,Ny)
    quad_f1 = zeros(Nx,Ny)
    quad_f2 = h./(3*ηav)

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
    @assert size(ηav)==(Nx,Ny)

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)


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