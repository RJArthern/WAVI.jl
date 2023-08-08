struct HGrid{T <: Real, N  <: Integer}
                  nxh :: N                                     # Number of grid cells in x-direction in HGrid
                  nyh :: N                                     # Number of grid cells in y-direction in HGrid
                 mask :: Array{Bool,2}                         # Mask specifying the model domain
            h_isfixed :: Array{Bool,2}                         # Mask specifying locations of fixed thickness
                    n :: N                                     # Total number of cells in the model domain
                 crop :: Diagonal{T,Array{T,1}}                # Crop matrix: diagonal matrix with mask entries on diag
                 samp :: SparseMatrixCSC{T,N}                  # Sampling matrix: take full domain to model domain 
               spread :: SparseMatrixCSC{T,N}                  # Sparse form of the sampling matrix 
              cent_xy :: KronType{T,N}                         # Centering operator from H-grid to C-grid
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
            nxh, 
            nyh,
            mask = trues(nxh,nyh),
            h_isfixed = falses(nxh,nxy),
            b,
            h,
            ηav = zeros(nxh,nyh),
            grounded_fraction = ones(nxh,nyh))

Construct a WAVI.jl HGrid with size (nxh,nyh)
HGrid stores fields that are defined on the problem's H grid. 
(Co-ordinates of HGrid stored in a Grid under xxh, yyh fields)

Keyword arguments
=================
    - 'nxh': (required) Number of grid cells in x-direction in HGrid (should be same as grid.nx)
            Note that we store the grid size here, even though it can be easily inferred from grid, to increase transparency in velocity solve.
    - 'nyh': (required) Number of grid cells in y-direction in HGrid (should be same as grid.ny)
    - 'mask': Mask specifying the model domain
    - 'h_isfixed': Mask specifying points where ice thickness is fixed
    - 'b': (requried) Bed elevation (bottom bathymetry)
    - 'h': (required) initial thickness of the ice
    - 'ηav': depth averaged visosity initially
    - 'grounded_fraction': initial grounded fraction
"""


function HGrid(;
                nxh, 
                nyh,
                mask = trues(nxh,nyh),
                h_isfixed = falses(nxh,nxy),
                b,
                h = zeros(nxh,nyh),
                ηav = zeros(nxh,nyh),
                grounded_fraction = ones(nxh,nyh))

    #check the sizes of inputs
    (size(mask) == size(h_isfixed) == size(b) == size(h) == size(ηav) == size(grounded_fraction) == (nxh,nyh)) || throw(DimensionMismatch("Sizes of inputs to HGrid must all be equal to nxh x nyh (i.e. $nxh x $nyh)"))

    #construct operators
    n = count(mask)
    crop = Diagonal(float(mask[:]))
    samp = sparse(1:n,(1:(nxh*nyh))[mask[:]],ones(n),n,nxh*nyh)
    spread = sparse(samp')
    cent_xy = c(nyh-1) ⊗ c(nxh-1)
    dneghηav = Ref(crop*Diagonal(zeros(nxh*nyh))*crop)
    dimplicit = Ref(crop*Diagonal(zeros(nxh*nyh))*crop)
     
    #construct quantities not passed
    s = zeros(nxh,nyh)
    dhdt = zeros(nxh,nyh) 
    accumulation = zeros(nxh,nyh)
    basal_melt = zeros(nxh,nyh)
    haf = zeros(nxh,nyh)
    dsdh = ones(nxh,nyh)
    shelf_strain_rate = zeros(nxh,nyh)
    av_speed = zeros(nxh,nyh) 
    u = zeros(nxh,nyh) 
    v = zeros(nxh,nyh)
    us = zeros(nxh,nyh) 
    vs = zeros(nxh,nyh)
    ub = zeros(nxh,nyh) 
    vb= zeros(nxh,nyh)
    bed_speed = zeros(nxh,nyh)
    weertman_c = zeros(nxh,nyh)
    β = zeros(nxh,nyh)
    βeff = zeros(nxh,nyh)
    τbed = zeros(nxh,nyh)
    quad_f1 = zeros(nxh,nyh)
    quad_f2 = zeros(nxh,nyh)
    quad_f2[mask] = h[mask]./(3*ηav[mask])

    #check sizes of everything
    @assert size(mask)==(nxh,nyh); #@assert mask == clip(mask)
    @assert size(h_isfixed)==(nxh,nyh); 
    @assert n == count(mask)
    @assert crop == Diagonal(float(mask[:]))
    @assert samp == sparse(1:n,(1:(nxh*nyh))[mask[:]],ones(n),n,nxh*nyh)
    @assert spread == sparse(samp')
    @assert size(cent_xy) == ((nxh-1)*(nyh-1),nxh*nyh)
    @assert size(b)==(nxh,nyh)
    @assert size(h)==(nxh,nyh)
    @assert size(s)==(nxh,nyh)
    @assert size(dhdt)==(nxh,nyh)
    @assert size(accumulation)==(nxh,nyh)
    @assert size(basal_melt)==(nxh,nyh)
    @assert size(haf)==(nxh,nyh)
    @assert size(grounded_fraction)==(nxh,nyh)
    @assert size(dsdh)==(nxh,nyh)
    @assert size(shelf_strain_rate)==(nxh,nyh)
    @assert size(u)==(nxh,nyh)
    @assert size(v)==(nxh,nyh)
    @assert size(av_speed)==(nxh,nyh) 
    @assert size(ub)==(nxh,nyh)
    @assert size(vb)==(nxh,nyh)
    @assert size(us)==(nxh,nyh)
    @assert size(vs)==(nxh,nyh)
    @assert size(bed_speed)==(nxh,nyh)
    @assert size(weertman_c)==(nxh,nyh)
    @assert size(β)==(nxh,nyh)
    @assert size(βeff)==(nxh,nyh)
    @assert size(τbed)==(nxh,nyh)
    @assert size(quad_f1)==(nxh,nyh)
    @assert size(quad_f2)==(nxh,nyh)
    @assert size(ηav)==(nxh,nyh)

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)
    h_isfixed = convert(Array{Bool,2}, h_isfixed)


return HGrid(
            nxh,
            nyh,
            mask,
            h_isfixed,
            n,
            crop,
            samp, 
            spread,
            cent_xy,
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
