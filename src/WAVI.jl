module WAVI

#Useful packages
using LinearAlgebra, SparseArrays, LinearMaps, Parameters,
      IterativeSolvers, Interpolations, BenchmarkTools, Reexport,
      NetCDF, JLD2, Setfield, MAT, ImageFiltering

#Import functions so they can be modified in this module.
import LinearAlgebra: ldiv!
import SparseArrays: spdiagm, spdiagm_internal, dimlub
import Setfield: @set

#This module will export these functions and types, allowing basic use of the model.
export
    #Structures
    Model, Params, TimesteppingParams, Grid, SolverParams, InitialConditions, OutputParams, Simulation,

    #Simulation controls
    update_state!, timestep!, run_simulation!,
#Abstract types
abstract type AbstractModel{T <: Real, N <: Integer} end
abstract type AbstractPreconditioner{T <: Real, N <: Integer} end

#Type alias, just for abreviation
const KronType{T,N} = LinearMaps.KroneckerMap{T,Tuple{LinearMaps.WrappedMap{T,SparseMatrixCSC{T,N}},
                        LinearMaps.WrappedMap{T,SparseMatrixCSC{T,N}}}} where {T <: Real, N <: Integer}

#Concrete types

#Struct to hold model parameters.
#Format: fieldname::Type = default_value.
#T & N are type parameters, usually real numbers (e.g. Float64) and integers (e.g. Int64) respectively.
@with_kw struct Params{T <: Real, N <: Integer}
nx::N = 80 
ny::N = 10
nσ::N = 4
x0::T = 0.0
y0::T = -40000.0
dx::T = 8000.0
dy::T = 8000.0
h_mask::Array{Bool,2}=trues(nx,ny);@assert size(h_mask)==(nx,ny);@assert h_mask == clip(h_mask)
u_iszero::Array{Bool,2}=falses(nx+1,ny);@assert size(u_iszero)==(nx+1,ny)
v_iszero::Array{Bool,2}=falses(nx,ny+1);@assert size(v_iszero)==(nx,ny+1)
dt::T = 1.0
g::T = 9.81
density_ice::T = 918.0
density_ocean::T = 1028.0
gas_const=8.314;
sec_per_year::T = 3600*24*365.25
starting_viscosity::T = 1.0e7
starting_temperature::T = 265.700709
bed_elevation::Array{T,2} = zeros(nx,ny); @assert size(bed_elevation)==(nx,ny)  
starting_thickness::Array{T,2} = zeros(nx,ny); @assert size(starting_thickness)==(nx,ny)
accumulation_rate::T = 0.0
basal_melt_rate::Array{T,2} = zeros(nx,ny); @assert size(basal_melt_rate)==(nx,ny)
starting_damage::T = 0.0
glen_a_activation_energy::T = 5.8631e+04
glen_a_ref::T = 4.9e-16 *sec_per_year * 1.0e-9
glen_temperature_ref::T = 263.15
glen_n::T = 3.0
glen_reg_strain_rate::T = 1.0e-5
n_iter_viscosity::N = 2;  @assert n_iter_viscosity ==2
weertman_c::T = 1e4    
weertman_m::T = 3.0
weertman_reg_speed::T = 1.0e-5
sea_level_wrt_geoid::T = 0.0
minimum_thickness::T = 50.0  
tol_picard::T = 1e-5
maxiter_picard::N = 30
tol_coarse::T = 1e-5
maxiter_coarse::N = 1000
levels::N = 3
wavelet_threshold::T = 10.0
nsmooth::N = 5
smoother_omega::T = 1.0
stencil_margin::N = 3
end

#Struct to hold information on h-grid, located at cell centers.
@with_kw struct HGrid{T <: Real, N <: Integer}
nx::N 
ny::N
x0::T = 0.0
y0::T = 0.0    
dx::T = 1.0
dy::T = 1.0
xx::Array{T,2}=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]; @assert size(xx)==(nx,ny)
yy::Array{T,2}=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]; @assert size(yy)==(nx,ny)
mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny); @assert mask == clip(mask)
n::N = count(mask); @assert n == count(mask)
crop::Diagonal{T,Array{T,1}} = Diagonal(float(mask[:])); @assert crop == Diagonal(float(mask[:]))
samp::SparseMatrixCSC{T,N} = crop[mask[:],:]; @assert samp == crop[mask[:],:]
spread::SparseMatrixCSC{T,N} = sparse(samp'); @assert spread == sparse(samp')
b::Array{T,2} = params.bed_elevation; @assert size(b)==(nx,ny)
h::Array{T,2} = zeros(nx,ny); @assert size(h)==(nx,ny)
s::Array{T,2} = zeros(nx,ny); @assert size(s)==(nx,ny)
dhdt::Array{T,2} = zeros(nx,ny); @assert size(dhdt)==(nx,ny) 
accumulation::Array{T,2} = zeros(nx,ny); @assert size(accumulation)==(nx,ny)
basal_melt::Array{T,2} = zeros(nx,ny); @assert size(basal_melt)==(nx,ny)   
haf::Array{T,2} = zeros(nx,ny); @assert size(haf)==(nx,ny)
grounded_fraction::Array{T,2} = ones(nx,ny); @assert size(grounded_fraction)==(nx,ny)
dsdh::Array{T,2} = ones(nx,ny); @assert size(dsdh)==(nx,ny)    
shelf_strain_rate::Array{T,2} = zeros(nx,ny); @assert size(shelf_strain_rate)==(nx,ny)
av_speed::Array{T,2} = zeros(nx,ny); @assert size(av_speed)==(nx,ny)
bed_speed::Array{T,2} = zeros(nx,ny); @assert size(bed_speed)==(nx,ny)
weertman_c::Array{T,2} = zeros(nx,ny); @assert size(weertman_c)==(nx,ny)
β::Array{T,2} = zeros(nx,ny); @assert size(β)==(nx,ny)
βeff::Array{T,2} = zeros(nx,ny); @assert size(βeff)==(nx,ny)
τbed::Array{T,2} = zeros(nx,ny); @assert size(τbed)==(nx,ny)
ηav::Array{T,2}; @assert size(ηav)==(nx,ny)
quad_f2::Array{T,2} = h./(3*ηav); @assert size(quad_f2)==(nx,ny)
dneghηav::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(zeros(nx*ny))*crop)
dimplicit::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(zeros(nx*ny))*crop)
end

#Struct to hold information on u-grid, located at grid-East and grid-West cell faces.
@with_kw struct UGrid{T <: Real, N <: Integer}
nx::N
ny::N
x0::T = 0.0
y0::T = 0.0    
dx::T = 1.0
dy::T = 1.0
xx::Array{T,2}=[x0+(i-1.0)*dx for i=1:nx, j=1:ny]; @assert size(xx)==(nx,ny)
yy::Array{T,2}=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]; @assert size(yy)==(nx,ny)
mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny)
n::N = count(mask); @assert n == count(mask)
crop::Diagonal{T,Array{T,1}} = Diagonal(float(mask[:])); @assert crop == Diagonal(float(mask[:]))
samp::SparseMatrixCSC{T,N} = crop[mask[:],:]; @assert samp == crop[mask[:],:]
spread::SparseMatrixCSC{T,N} = sparse(samp'); @assert spread == sparse(samp')
cent::KronType{T,N} = spI(ny) ⊗ c(nx-1)
∂x::KronType{T,N} = spI(ny) ⊗ ∂1d(nx-1,dx)
∂y::KronType{T,N} = ∂1d(ny-1,dy) ⊗ χ(nx-2)
levels::N     
dwt::KronType{T,N} = wavelet_matrix(ny,levels,"forward" ) ⊗ wavelet_matrix(nx,levels,"forward")
s::Array{T,2} = zeros(nx,ny); @assert size(s)==(nx,ny)
h::Array{T,2} = zeros(nx,ny); @assert size(h)==(nx,ny)
grounded_fraction::Array{T,2} = ones(nx,ny); @assert size(grounded_fraction)==(nx,ny)
βeff::Array{T,2} = zeros(nx,ny); @assert size(βeff)==(nx,ny)
dnegβeff::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(-βeff[:])*crop)
u::Array{T,2} = zeros(nx,ny); @assert size(u)==(nx,ny)
end

#Struct to hold information on v-grid, located at grid-North and grid-South cell faces.
@with_kw struct VGrid{T <: Real, N <: Integer}
nx::N
ny::N
x0::T = 0.0
y0::T = 0.0    
dx::T = 1.0
dy::T = 1.0
xx::Array{T,2}=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]; @assert size(xx)==(nx,ny)
yy::Array{T,2}=[y0+(j-1.0)*dy for i=1:nx, j=1:ny]; @assert size(yy)==(nx,ny)
mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny)
n::N = count(mask); @assert n == count(mask)
crop::Diagonal{T,Array{T,1}} = Diagonal(float(mask[:])); @assert crop == Diagonal(float(mask[:]))
samp::SparseMatrixCSC{T,N} = crop[mask[:],:]; @assert samp == crop[mask[:],:]
spread::SparseMatrixCSC{T,N} = sparse(samp'); @assert spread == sparse(samp')
cent::KronType{T,N} = c(ny-1) ⊗ spI(nx)
∂x::KronType{T,N} = χ(ny-2) ⊗ ∂1d(nx-1,dx)
∂y::KronType{T,N} = ∂1d(ny-1,dy) ⊗ spI(nx)
levels::N  
dwt::KronType{T,N} = wavelet_matrix(ny,levels,"forward" ) ⊗ wavelet_matrix(nx,levels,"forward")
s::Array{T,2} = zeros(nx,ny); @assert size(s)==(nx,ny)
h::Array{T,2} = zeros(nx,ny); @assert size(h)==(nx,ny)
grounded_fraction::Array{T,2} = ones(nx,ny); @assert size(grounded_fraction)==(nx,ny)
βeff::Array{T,2} = zeros(nx,ny); @assert size(βeff)==(nx,ny)
dnegβeff::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(crop*Diagonal(-βeff[:])*crop)
v::Array{T,2} = zeros(nx,ny); @assert size(v)==(nx,ny)
end

#Struct to hold information on c-grid, located at cell corners.
@with_kw struct CGrid{T <: Real, N <: Integer}
nx::N
ny::N
x0::T = 0.0
y0::T = 0.0    
dx::T = 1.0
dy::T = 1.0
xx::Array{T,2}=[x0+i*dx for i=1:nx, j=1:ny]; @assert size(xx)==(nx,ny)
yy::Array{T,2}=[y0+j*dy for i=1:nx, j=1:ny]; @assert size(yy)==(nx,ny)
mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny)
n::N = count(mask); @assert n == count(mask)
crop::Diagonal{T,Array{T,1}} = Diagonal(float(mask[:])); @assert crop == Diagonal(float(mask[:]))
samp::SparseMatrixCSC{T,N} = crop[mask[:],:]; @assert samp == crop[mask[:],:]
spread::SparseMatrixCSC{T,N} = sparse(samp'); @assert spread == sparse(samp')
cent::KronType{T,N} = sparse(c(ny)') ⊗ sparse(c(nx)')
end

#Struct to hold information on 3d-grid, (extends h-grid to multiple sigma levels).
@with_kw struct SigmaGrid{T <: Real, N <: Integer}
nx::N
ny::N
nσ::N
σ::Vector{T} = collect(range(0.0,length=nσ,stop=1.0)); @assert length(σ) == nσ
ζ::Vector{T} = one(eltype(σ)) .- σ ; @assert length(ζ) == nσ
quadrature_weights::Vector{T} = [0.5;ones(nσ-2);0.5]/(nσ-1); @assert length(quadrature_weights) == nσ
η::Array{T,3} = fill(params.starting_viscosity,nx,ny,nσ); @assert size(η)==(nx,ny,nσ)
θ::Array{T,3} = fill(params.starting_temperature,nx,ny,nσ); @assert size(θ)==(nx,ny,nσ)
Φ::Array{T,3} = fill(params.starting_damage,nx,ny,nσ); @assert size(Φ)==(nx,ny,nσ)
glen_b::Array{T,3} = glen_b.(θ,Φ); @assert size(glen_b)==(nx,ny,nσ)    
end

#Struct to hold information on wavelet-grid (u-component). 
@with_kw struct UWavelets{T <: Real, N <: Integer}
nx::N
ny::N
mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny)
n::Base.RefValue{N} = Ref(count(mask)); @assert n[] == count(mask)
crop::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(Diagonal(float(mask[:]))); @assert crop[] == Diagonal(float(mask[:]))
samp::Base.RefValue{SparseMatrixCSC{T,N}} = Ref(crop[][mask[:],:]); @assert samp[] == crop[][mask[:],:]
spread::Base.RefValue{SparseMatrixCSC{T,N}} = Ref(sparse(samp[]')); @assert spread[] == sparse(samp[]')
levels::N
idwt::KronType{T,N} = wavelet_matrix(ny,levels,"reverse" ) ⊗ wavelet_matrix(nx,levels,"reverse")
wavelets::Array{T,2} = zeros(nx,ny); @assert size(wavelets)==(nx,ny)
end

#Struct to hold information on wavelet-grid (v-component). 
@with_kw struct VWavelets{T <: Real, N <: Integer}
nx::N
ny::N
mask::Array{Bool,2} = trues(nx,ny); @assert size(mask)==(nx,ny)
n::Base.RefValue{N} = Ref(count(mask)); @assert n[] == count(mask)
crop::Base.RefValue{Diagonal{T,Array{T,1}}} = Ref(Diagonal(float(mask[:]))); @assert crop[] == Diagonal(float(mask[:]))
samp::Base.RefValue{SparseMatrixCSC{T,N}} = Ref(crop[][mask[:],:]); @assert samp[] == crop[][mask[:],:]
spread::Base.RefValue{SparseMatrixCSC{T,N}} = Ref(sparse(samp[]')); @assert spread[] == sparse(samp[]')
levels::N
idwt::KronType{T,N} = wavelet_matrix(ny,levels,"reverse" ) ⊗ wavelet_matrix(nx,levels,"reverse")
wavelets::Array{T,2} = zeros(nx,ny); @assert size(wavelets)==(nx,ny)
end

#Struct to hold model state comprised of all the above information.
@with_kw struct State{T <: Real, N <: Integer} <: AbstractModel{T,N}
params::Params{T,N}
gh::HGrid{T,N}
gu::UGrid{T,N}
gv::VGrid{T,N}
gc::CGrid{T,N}
g3d::SigmaGrid{T,N}   
wu::UWavelets{T,N}
wv::VWavelets{T,N}
end

#Struct to hold information about wavelet-based multigrid preconditioner.
@with_kw struct Preconditioner{T <: Real, N <: Integer} <: AbstractPreconditioner{T,N}
op::LinearMap{T}
op_diag::Vector{T} = diag(sparse(op))
nsmooth::N = 5    
sweep::Vector{N}
sweep_order::Vector{N} = unique(sweep)
smoother_omega::T = 1.0
restrict::LinearMap{T}
prolong::LinearMap{T}
op_coarse::LinearMap{T} = restrict*op*prolong
correction_coarse::Vector{T} = zeros(T,size(op_coarse,2))
tol_coarse::T = 1e-7
maxiter_coarse::N = 1000    
end


#Functions


#1D Matrix operator utility functions.
spI(n) = spdiagm(n,n, 0 => ones(n))
∂1d(n,dx) = spdiagm(n,n+1,0 => -ones(n), 1 => ones(n))/dx
c(n) = spdiagm(n,n+1,0 => ones(n), 1 => ones(n))/2
χ(n) = spdiagm(n,n+2, 1 => ones(n))

"""
    start(params)

Create WAVI State from input parameters.
"""
function start(params)
      
    #Define masks for points on h-, u-, v- and c-grids that lie in model domain.
    @assert params.h_mask==clip(params.h_mask) "Model domain mask has invalid points. Use clip function to remove them."
    h_mask = params.h_mask
    u_mask = get_u_mask(h_mask)
    v_mask = get_v_mask(h_mask)
    c_mask = get_c_mask(h_mask)
    
    #Remove all points on u- and v-grids with homogenous Dirichlet conditions.
    u_mask[params.u_iszero].=false
    v_mask[params.v_iszero].=false
      
    #h-grid
    gh=HGrid(x0=params.x0,
    y0=params.y0,        
    nx=params.nx,
    ny=params.ny,
    dx=params.dx,
    dy=params.dy,
    mask=h_mask,
    b = params.bed_elevation,
    h = params.starting_thickness,
    ηav = fill(params.starting_viscosity,params.nx,params.ny),
    )    
    
    #u-grid
    gu=UGrid(
    x0=params.x0,
    y0=params.y0,        
    nx=params.nx+1,
    ny=params.ny,
    dx=params.dx,
    dy=params.dy,
    mask=u_mask,
    levels=params.levels
    )  
    
    #v-grid
    gv=VGrid(
    x0=params.x0,
    y0=params.y0,        
    nx=params.nx,
    ny=params.ny+1,
    dx=params.dx,
    dy=params.dy,
    mask=v_mask,
    levels=params.levels
    ) 
    
    #c-grid
    gc=CGrid(
    x0=params.x0,
    y0=params.y0,        
    nx=params.nx-1,
    ny=params.ny-1,
    dx=params.dx,
    dy=params.dy,
    mask=c_mask
    ) 
    
    #3D-grid
    g3d=SigmaGrid(
    nx=params.nx,
    ny=params.ny,
    nσ=params.nσ,
    η = fill(params.starting_viscosity,params.nx,params.ny,params.nσ),
    θ = fill(params.starting_temperature,params.nx,params.ny,params.nσ),
    Φ = fill(params.starting_damage,params.nx,params.ny,params.nσ),
    glen_b = fill(glen_b(params.starting_temperature,params.starting_damage,params),params.nx,params.ny,params.nσ)    
    )
    
    #Wavelet-grid, u-component.
    wu=UWavelets(nx=params.nx+1,ny=params.ny,levels=params.levels)
    
    #Wavelet-grid, v-component.
    wv=VWavelets(nx=params.nx,ny=params.ny+1,levels=params.levels)
    
    #Use type constructor to build initial state.
    wavi=State(params,gh,gu,gv,gc,g3d,wu,wv)

    return wavi
end    

"""
    run!(wavi)

Run WAVI model for one timestep. Input variable is modified.
"""
function run!(wavi)
    update_surface_elevation!(wavi)
    update_geometry_on_uv_grids!(wavi)
    update_height_above_floatation!(wavi)
    update_grounded_fraction_on_huv_grids!(wavi)
    update_accumulation_rate!(wavi)
    update_basal_melt!(wavi)
    update_weertman_c!(wavi)
    update_dsdh!(wavi)
    update_velocities!(wavi)
    update_dhdt!(wavi)
    update_thickness!(wavi)
    update_wavelets!(wavi)
    return wavi
end

# Utility functions

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
    glen_b(temperature,damage,params)

Compute stiffness parameter B in Glen flow law.

"""

function glen_b(temperature,damage,params)
    @unpack glen_a_ref, glen_n, glen_a_activation_energy, glen_temperature_ref, gas_const = params
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



# Functions to update various quantities

"""
    update_geometry_on_uv_grids!(wavi::AbstractModel)

Interpolate thickness and surface elvation from h-grid to u- and v-grids. 

"""
function update_geometry_on_uv_grids!(wavi::AbstractModel)
    @unpack gh,gu,gv,gc=wavi
    onesvec=ones(gh.nx*gh.ny)
    gu.h[gu.mask].=(gu.samp*(gu.cent'*(gh.crop*gh.h[:])))./(gu.samp*(gu.cent'*(gh.crop*onesvec)))
    gu.s[gu.mask].=(gu.samp*(gu.cent'*(gh.crop*gh.s[:])))./(gu.samp*(gu.cent'*(gh.crop*onesvec)))
    gv.h[gv.mask].=(gv.samp*(gv.cent'*(gh.crop*gh.h[:])))./(gv.samp*(gv.cent'*(gh.crop*onesvec)))
    gv.s[gv.mask].=(gv.samp*(gv.cent'*(gh.crop*gh.s[:])))./(gv.samp*(gv.cent'*(gh.crop*onesvec)))
    return wavi
end

"""
    update_velocities!(wavi::AbstractModel)

 Solve momentum equation to update the velocities, plus Picard iteration for non-linear rheology.

"""
function update_velocities!(wavi::AbstractModel)
    @unpack gu,gv,wu,wv,params=wavi
    
    n = gu.n + gv.n  
    
    x=get_start_guess(wavi)    
    b=get_rhs(wavi)        
    
    rel_resid=zero(eltype(b))
    
    converged::Bool = false
    i_picard::Int64 = 0
    while !converged && (i_picard < params.maxiter_picard)
        
        i_picard = i_picard + 1
        
        set_velocities!(wavi,x)
        update_shelf_strain_rate!(wavi)
        update_av_speed!(wavi)
        update_bed_speed!(wavi)
        update_β!(wavi)
        update_basal_drag!(wavi)
        inner_update_viscosity!(wavi)
        update_av_viscosity!(wavi)
        update_quadrature_f2!(wavi)
        update_βeff!(wavi)
        update_βeff_on_uv_grids!(wavi)
        update_rheological_operators!(wavi)
        op=get_op(wavi)
        
        rel_resid = norm(b .- op*x)/norm(b)                
        converged = rel_resid < params.tol_picard
         
        p=get_preconditioner(wavi,op)   
        precondition!(x, p, b)
        
    end  
    set_velocities!(wavi,x)
    
    return wavi
end
"""
    get_start_guess(wavi::AbstractModel)

 Return starting guess used to begin iterative solution of velocities.

"""
function get_start_guess(wavi::AbstractModel)
    @unpack gu,gv=wavi
    @assert eltype(gu.u)==eltype(gv.v)
    n = gu.n + gv.n
    x=[gu.samp*gu.u[:];gv.samp*gv.v[:]]
    return x
end

"""
    get_rhs(wavi::AbstractModel)

 Return right hand side vector of momentum equations.

"""

function get_rhs(wavi::AbstractModel)
    @unpack gh,gu,gv,gc,params=wavi
    onesvec=ones(gh.nx*gh.ny)
    surf_elev_adjusted = gh.crop*(gh.s[:] .+ params.dt*gh.dsdh[:].*(gh.accumulation[:].-gh.basal_melt[:]))
    f1=[
        (params.density_ice*params.g*gu.h[gu.mask]).*(gu.samp*(-gu.∂x'*surf_elev_adjusted))
        ;
        (params.density_ice*params.g*gv.h[gv.mask]).*(gv.samp*(-gv.∂y'*surf_elev_adjusted))
       ]
    f2=[
        (0.5*params.density_ice*params.g*gu.h[gu.mask].^2
        -0.5*params.density_ocean*params.g*(icedraft.(gu.s[gu.mask],gu.h[gu.mask],params.sea_level_wrt_geoid)).^2
        -params.density_ice*params.g*gu.h[gu.mask].*gu.s[gu.mask]).*gu.samp*(-gu.∂x'*(gh.crop*onesvec))
        ;
        (0.5*params.density_ice*params.g*gv.h[gv.mask].^2
        -0.5*params.density_ocean*params.g*(icedraft.(gv.s[gv.mask],gv.h[gv.mask],params.sea_level_wrt_geoid)).^2
        -params.density_ice*params.g*gv.h[gv.mask].*gv.s[gv.mask]).*gv.samp*(-gv.∂y'*(gh.crop*onesvec))
        ]
    rhs=f1+f2
    return rhs
end

"""
    get_op(wavi::AbstractModel{T,N}) where {T,N}

 Get operator, defined as a LinearMap type.

"""
function get_op(wavi::AbstractModel{T,N}) where {T,N}
    @unpack gu,gv=wavi
    n = gu.n + gv.n
    op_fun(x) = opvec(wavi,x)
    op=LinearMap{T}(op_fun,n;issymmetric=true,ismutating=false,ishermitian=true,isposdef=true)
end
"""
    get_op_diag(wavi::AbstractModel,op::LinearMap)

 Get diagonal of operator for use in preconditioner.

"""
function get_op_diag(wavi::AbstractModel,op::LinearMap)
    @unpack gu,gv,params=wavi
    m,n=size(op)
    @assert m == n == gu.n + gv.n
    op_diag=zeros(eltype(op),n)
    sm=params.stencil_margin
    sweep=[[1+mod((i-1),sm)+sm*mod((j-1),sm) for i=1:gu.nx, j=1:gu.ny][gu.mask];
           [1+sm^2+mod((i-1),sm)+sm*mod((j-1),sm) for i=1:gv.nx, j=1:gv.ny][gv.mask] ]
    e=zeros(Bool,n)
    for i = unique(sweep)
        e[sweep .== i] .= true
        op_diag[e] .= (op*e)[e]
        e[sweep .== i] .= false
    end
    return op_diag
end
"""
    get_preconditioner(wavi::AbstractModel{T,N},op::LinearMap{T}) where {T, N}

 Get precondtioner. Details are stored in a struct of type Preconditioner that can be passed to a solver.

"""
function get_preconditioner(wavi::AbstractModel{T,N},op::LinearMap{T}) where {T, N}

    @unpack gu,gv,wu,wv,params=wavi
    
    m,n=size(op)
    @assert m == n == gu.n + gv.n
    
    n = gu.n + gv.n
    n_coarse = wu.n[] + wv.n[]
    
    restrict_fun(x) = restrictvec(wavi,x)
    prolong_fun(x) = prolongvec(wavi,x)
    
    restrict=LinearMap{T}(restrict_fun,n_coarse,n;issymmetric=false,ismutating=false,ishermitian=false,isposdef=false)
    prolong=LinearMap{T}(prolong_fun,n,n_coarse;issymmetric=false,ismutating=false,ishermitian=false,isposdef=false)    
    
    op_diag=get_op_diag(wavi,op)
    
    #Four colour Jacobi preconditioner. Red-Black checkerboard Jacobi for each velocity component.
    sweep=[[1 .+ mod(i-j,2) for i=1:gu.nx, j=1:gu.ny][gu.mask];[3 .+ mod(i-j,2) for i=1:gv.nx, j=1:gv.ny][gv.mask]]
    sweep_order=[1,3,2,4]
    
    p=Preconditioner{T,N}(op=op, restrict=restrict, prolong=prolong,sweep=sweep, sweep_order=sweep_order, 
            op_diag=op_diag, nsmooth=params.nsmooth, tol_coarse = params.tol_coarse, 
            maxiter_coarse = params.maxiter_coarse, smoother_omega=params.smoother_omega)
    
    return p
end
"""
    LinearAlgebra.ldiv!(x::AbstractVecOrMat{T}, p::AbstractPreconditioner{T,N}, b::AbstractVecOrMat{T}) where {T,N}

    #Melt ratesS
    PlumeEmulator, BinfileMeltRate, UniformMeltRate, MISMIPMeltRateOne, PICO, QuadraticMeltRate, UniformMeltFloatOnly, UniformMeltFloatOnlyBasinSpecific, QuadraticMeltRateExponentVariation, QuadraticMeltRateExponentVariationBasins, MeltRateExponentVariationBasins, 

    #Post-processing controls
    volume_above_floatation, height_above_floatation

#Reexport Modules useful for users of the WAVI module
@reexport using JLD2
@reexport using Setfield

#Abstract types
abstract type AbstractGrid{T <: Real, N <: Integer} end
abstract type AbstractMeltRate end
abstract type AbstractModel{T <: Real, N <: Integer, M <: AbstractMeltRate} end
abstract type AbstractPreconditioner{T <: Real, N <: Integer} end
#abstract type AbstractSimulation{T,N,R,A,W} end



#Type alias, just for abreviation
const KronType{T,N} = LinearMaps.KroneckerMap{T,Tuple{LinearMaps.WrappedMap{T,SparseMatrixCSC{T,N}},
                        LinearMaps.WrappedMap{T,SparseMatrixCSC{T,N}}}} where {T <: Real, N <: Integer}

#Concrete types

##################################################################################
#include all of the code
include("OutputParams/OutputParams.jl")
include("Grid.jl")
include("Params.jl")
include("SolverParams.jl")
include("TimesteppingParams.jl")
include("Clock.jl")
include("InitialConditions.jl")
include("Wavelets/Wavelets.jl")
include("Fields/Fields.jl")
include("Models/Model.jl")
include("MeltRate/MeltRate.jl")
include("Simulations/Simulation.jl")
include("utilities.jl")


"""
    spdiagm(m::Integer, n::Integer, kv::Pair{<:Integer,<:AbstractVector}...)

Method from Julia V1.4.2 to create non-square sparse matrix from diagonals. Included for backward compatibility.
"""
spdiagm(m::Integer, n::Integer, kv::Pair{<:Integer,<:AbstractVector}...) = _spdiagm((Int(m),Int(n)), kv...)
function _spdiagm(size, kv::Pair{<:Integer,<:AbstractVector}...)
    I, J, V = spdiagm_internal(kv...)
    mmax, nmax = dimlub(I), dimlub(J)
    mnmax = max(mmax, nmax)
    m, n = something(size, (mnmax,mnmax))
    (m ≥ mmax && n ≥ nmax) || throw(DimensionMismatch("invalid size=$size"))
    return sparse(I, J, V, m, n)
end
end


