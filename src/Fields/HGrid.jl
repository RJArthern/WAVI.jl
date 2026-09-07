struct HGrid{T <: Real, N  <: Integer}
                  nxh :: N                                     # Number of grid cells in x-direction in HGrid
                  nyh :: N                                     # Number of grid cells in y-direction in HGrid
                 mask :: Array{Bool,2}                         # Mask specifying the model domain
            h_isfixed :: Array{Bool,2}                         # Mask specifying locations of fixed thickness
hyd_potential_isfixed :: Array{Bool,2}                         # Mask specifying locations of fixed hydraulic potential at the bed
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
           basal_melt :: Array{T,2}                            # basal melt rate  
                  haf :: Array{T,2}                            # Grid cell height above floatation
    grounded_fraction :: Array{T,2}                            # Grid cell grounded fraction 
                 dsdh :: Array{T,2}                            # Change of surface elevation per unit thickness change
    shelf_strain_rate :: Array{T,2}                            # Strain rate appropriate for shelf (no basal drag) 
             av_speed :: Array{T,2}                            # Depth averaged speed 
                    u :: Array{T,2}                            # Depth averaged x-velocity 
                    v :: Array{T,2}                            # Depth averaged y-velocity 
                   us :: Array{T,2}                            # x-velocity at the surface 
                   vs :: Array{T,2}                            # y-velocity at the surface 
                   ub :: Array{T,2}                            # x-velocity at the bed 
                   vb :: Array{T,2}                            # y-velocity at the bed
            bed_speed :: Array{T,2}                            # Ice speed at the bed
     drag_coefficient :: Array{T,2}                            # Sliding law drag coefficients
                    β :: Array{T,2}                            # Raw β value (eqn 8 in Arthern 2015 JGeophysRes)
                 βeff :: Array{T,2}                            # Effective β value (eqn 12 in Arthern 2015 JGeophysRes)
                 τbed :: Array{T,2}                            # Stress at the bed
                  ηav :: Array{T,2}                            # Depth averaged viscosity
              quad_f0 :: Array{T,2}                            # F0 quadratrue field (eqn 7 in Arthern 2015 JGeophysRes)
              quad_f1 :: Array{T,2}                            # F1 quadratrue field (eqn 7 in Arthern 2015 JGeophysRes)
              quad_f2 :: Array{T,2}                            # F2 quadrature field (eqn 7 in Arthern 2015 JGeophysRes)
             dneghηav :: Base.RefValue{Diagonal{T,Array{T,1}}} # Rheological operator (-h × ηav)
            dimplicit :: Base.RefValue{Diagonal{T,Array{T,1}}} # Rheological operator (-ρi × g × dt × dshs)
basal_water_thickness :: Array{T,2}                            # basal water thickness
hydraulic_potential_b :: Array{T,2}                            # hydraulic potential at the bed
   effective_pressure :: Array{T,2}                            # effective pressure  
     shelf_basal_melt :: Array{T,2}                            # basal melt rate under ice shelves (ie floating ice)
                θ_ave :: Array{T,2}                            # depth-averaged temperature
              σzzsurf :: Array{T,2}                            # Sigmazzsurf calculated in inversion (see in Arthern 2015 JGeophysRes)
              τx_surf :: Array{T,2}                            # Stress at surface u component
              τy_surf :: Array{T,2}                            # Stress at surface u component
              τsurf :: Array{T,2}                              # Stress at surface
              τx_bed :: Array{T,2}                             # Stress at bed u component
              τy_bed :: Array{T,2}                             # Stress at bed v component
           surf_speed :: Array{T,2}                            # Ice speed at the surface u magnitude
        shelf_heating :: Array{T,2}                            # Shelf heating rate calculated in inversion  (see in Arthern 2015 JGeophysRes)
   vert_shear_heating :: Array{T,2}                            # Vertical shear heating rate calculated in inversion (see in Arthern 2015 JGeophysRes)
         drag_heating :: Array{T,2}                            # SDrag heating rate calculated in inversion (see in Arthern 2015 JGeophysRes)
         preBfactor   :: Array{T,2}                            # preBfactor calculated in inversion (see in Arthern 2015 JGeophysRes)
             mpi_rank :: Array{T,2}                            # MPI rank for each cell (for visualization)
end


"""
    HGrid(;
            nxh, 
            nyh,
            mask = trues(nxh,nyh),
            h_isfixed = falses(nxh,nyh),
            hyd_potential_isfixed = falses(nxh,nyh),
            b,
            h,
            ηav = zeros(nxh,nyh),
            grounded_fraction = ones(nxh,nyh),
            basal_water_thickness = zeros(nxh,nyh),
            hydraulic_potential_b = zeros(nxh,nyh)
	        effective_pressure = zeros(nxh,nyh),
            basal_melt = zeros(nxh,nyh),
            θ_ave = zeros(nxh,nyh),
            grounded_fraction = ones(nxh,nyh),
            preBfactor = ones(nxh,nyh),
            mpi_rank = -ones(nxh,nyh))

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
    - 'hyd_potential_isfixed': Mask specifying points where the hydraulic potential at the bed is fixed
    - 'b': (required) Bed elevation (bottom bathymetry)
    - 'h': initial thickness of the ice
    - 'ηav': depth averaged visosity initially
    - 'grounded_fraction': initial grounded fraction
    - 'basal_water_thickness' : initial basal water thickness
    - 'hydraulic_potential_b' : initial hydraulic potential at the bed
    - 'effective_pressure': initial effective pressure
    - 'basal_melt': initial basal melt rate
    - 'θ_ave': initial depth-averaged temperature
    - 'preBfactor : preBfactor (1=no viscosity enhancement)
    - 'mpi_rank: field to store MPI rank identifier  
"""


function HGrid(;
                nxh, 
                nyh,
                mask = trues(nxh,nyh),
                h_isfixed = falses(nxh,nyh),
                hyd_potential_isfixed = falses(nxh,nyh),
                b,
                h = zeros(nxh,nyh),
                ηav = zeros(nxh,nyh),
                grounded_fraction = ones(nxh,nyh),
                basal_water_thickness = zeros(nxh,nyh),
                hydraulic_potential_b = zeros(nxh,nyh),
		        effective_pressure = zeros(nxh,nyh),
                basal_melt = zeros(nxh,nyh),
                θ_ave = zeros(nxh,nyh),
                preBfactor = ones(nxh,nyh),
                mpi_rank = -ones(nxh,nyh),
                storage_only::Bool = false,
)

    @debug "HGrid input sizes" sizes = (
        size_mask = size(mask),
        size_h_isfixed = size(h_isfixed),
        size_hyd_potential_isfixed = size(hyd_potential_isfixed),
        size_b = size(b),
        size_h = size(h),
        size_ηav = size(ηav),
        size_grounded_fraction = size(grounded_fraction),
        size_basal_water_thickness = size(basal_water_thickness),
        size_hydraulic_potential_b = size(hydraulic_potential_b),
        size_effective_pressure = size(effective_pressure),
        size_basal_melt = size(basal_melt),
        size_θ_ave = size(θ_ave),
        size_preBfactor = size(preBfactor),
        size_mpi_rank = size(mpi_rank),
    )

    #check the sizes of inputs
    (size(mask) == size(h_isfixed) == size(hyd_potential_isfixed) == size(b) == size(h) == size(ηav) == size(grounded_fraction) == size(basal_water_thickness) == size(hydraulic_potential_b) == size(effective_pressure) == size(basal_melt) == size(θ_ave) == size(preBfactor) == size(mpi_rank) == (nxh,nyh)) || throw(DimensionMismatch("Sizes of inputs to HGrid must all be equal to nxh x nyh (i.e. $nxh x $nyh)"))

    #construct operators
    if storage_only
        n = zero(nxh)
        crop = Diagonal(Float64[])
        samp = spzeros(Float64, 0, nxh * nyh)
        spread = spzeros(Float64, nxh * nyh, 0)
        cent_xy = _storage_only_kron()
        dneghηav = Ref(Diagonal(Float64[]))
        dimplicit = Ref(Diagonal(Float64[]))
    else
        n = count(mask)
        crop = Diagonal(float(mask[:]))
        samp = sparse(1:n,(1:(nxh*nyh))[mask[:]],ones(n),n,nxh*nyh)
        spread = sparse(samp')
        cent_xy = c(nyh-1) ⊗ c(nxh-1)
        dneghηav = Ref(crop*Diagonal(zeros(nxh*nyh))*crop)
        dimplicit = Ref(crop*Diagonal(zeros(nxh*nyh))*crop)
    end
     
    #construct quantities not passed
    s = zeros(nxh,nyh)
    dhdt = zeros(nxh,nyh) 
    accumulation = zeros(nxh,nyh)
    shelf_basal_melt = zeros(nxh,nyh)
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
    drag_coefficient = zeros(nxh,nyh)
    β = zeros(nxh,nyh)
    βeff = zeros(nxh,nyh)
    τbed = zeros(nxh,nyh)
    quad_f0 = zeros(nxh,nyh)
    quad_f1 = zeros(nxh,nyh)
    quad_f2 = zeros(nxh,nyh)
    if !storage_only
        quad_f2[mask] = h[mask]./(3*ηav[mask])
    end
    σzzsurf=zeros(nxh,nyh) 
    τx_surf=zeros(nxh,nyh) 
    τy_surf=zeros(nxh,nyh) 
    τsurf=zeros(nxh,nyh) 
    τx_bed=zeros(nxh,nyh) 
    τy_bed =zeros(nxh,nyh) 
    surf_speed=zeros(nxh,nyh) 
    shelf_heating=zeros(nxh,nyh) 
    vert_shear_heating=zeros(nxh,nyh) 
    drag_heating=zeros(nxh,nyh) 
    #
    #preBfactor=ones(nxh,nyh) 


    # mpi_rank is passed as a parameter (defaults to -1 for non-MPI runs)

    #check sizes of everything
    @assert size(mask)==(nxh,nyh); #@assert mask == clip(mask)
    @assert size(h_isfixed)==(nxh,nyh);
    @assert size(hyd_potential_isfixed)==(nxh,nyh);
    if !storage_only
        @assert n == count(mask)
        @assert size(cent_xy) == ((nxh-1)*(nyh-1),nxh*nyh)
    end
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
    @assert size(drag_coefficient)==(nxh,nyh)
    @assert size(β)==(nxh,nyh)
    @assert size(βeff)==(nxh,nyh)
    @assert size(τbed)==(nxh,nyh)
    @assert size(quad_f0)==(nxh,nyh)
    @assert size(quad_f1)==(nxh,nyh)
    @assert size(quad_f2)==(nxh,nyh)
    @assert size(ηav)==(nxh,nyh)
    @assert size(basal_water_thickness)==(nxh,nyh)
    @assert size(hydraulic_potential_b)==(nxh,nyh)
    @assert size(effective_pressure)==(nxh,nyh)
    @assert size(shelf_basal_melt)==(nxh,nyh)
    @assert size(θ_ave)==(nxh,nyh)
    @assert size(σzzsurf)==(nxh,nyh)
    @assert size(τx_surf)==(nxh,nyh)
    @assert size(τy_surf)==(nxh,nyh)
    @assert size(τsurf)==(nxh,nyh)
    @assert size(τx_bed)==(nxh,nyh)
    @assert size(τy_bed)==(nxh,nyh)
    @assert size(surf_speed)==(nxh,nyh)
    @assert size(shelf_heating)==(nxh,nyh)
    @assert size(vert_shear_heating)==(nxh,nyh)
    @assert size(drag_heating)==(nxh,nyh)
    @assert size(preBfactor)==(nxh,nyh)
    @assert size(mpi_rank)==(nxh,nyh)

    #make sure boolean type rather than bitarray
    mask = convert(Array{Bool,2}, mask)
    h_isfixed = convert(Array{Bool,2}, h_isfixed)
    hyd_potential_isfixed = convert(Array{Bool,2}, hyd_potential_isfixed)


return HGrid(
            nxh,
            nyh,
            mask,
            h_isfixed,
            hyd_potential_isfixed,
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
            drag_coefficient,
            β,
            βeff,
            τbed,
            ηav,
            quad_f0,
            quad_f1,
            quad_f2,
            dneghηav,
            dimplicit,
            basal_water_thickness,
            hydraulic_potential_b,
	        effective_pressure,
            shelf_basal_melt,
            θ_ave,
            σzzsurf, 
            τx_surf,
            τy_surf,
            τsurf,
            τx_bed,
            τy_bed, 
            surf_speed,
            shelf_heating,
            vert_shear_heating,
            drag_heating,
            preBfactor,
            mpi_rank
        )
end
