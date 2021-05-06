module WAVI

#Useful packages
using LinearAlgebra, SparseArrays, LinearMaps, Parameters,
      IterativeSolvers, Interpolations, BenchmarkTools, Reexport, NetCDF, JLD2, HDF5, Setfield

#Import functions so they can be modified in this module.
import LinearAlgebra: ldiv!
import SparseArrays: spdiagm, spdiagm_internal, dimlub
import Setfield: @set

#This module will export these functions and types, allowing basic use of the model.
export start, run!, plot_output, State, Params, 
TimesteppingParams, Grid, SolverParams, InitialConditions, simulation

#Reexport Modules useful for users of the WAVI module
@reexport using JLD2

#Abstract types
abstract type AbstractGrid{T <: Real, N <: Integer} end
abstract type AbstractModel{T <: Real, N <: Integer} end
abstract type AbstractPreconditioner{T <: Real, N <: Integer} end

#Type alias, just for abreviation
const KronType{T,N} = LinearMaps.KroneckerMap{T,Tuple{LinearMaps.WrappedMap{T,SparseMatrixCSC{T,N}},
                        LinearMaps.WrappedMap{T,SparseMatrixCSC{T,N}}}} where {T <: Real, N <: Integer}

#Concrete types



##################################################################################
#grid info
struct Grid{T <: Real, N <: Integer} <: AbstractGrid{T,N}
    nx::N
    ny::N
    nσ::N 
    dx::T
    dy::T
    x0::T
    y0::T
    h_mask::Array{Bool,2}
    u_iszero::Array{Bool,2} #zero boundary condition locations on u
    v_iszero::Array{Bool,2} #zero boundary condition locations on u
    xxh::Array{T,2}         #x co-ordinates matrix of h grid
    yyh::Array{T,2}         #y co-ordinates matrix of h grid
    xxu::Array{T,2}         #x co-ordinates matrix of u grid
    yyu::Array{T,2}         #y co-ordinates matrix of u grid
    xxv::Array{T,2}         #x co-ordinates matrix of v grid
    yyv::Array{T,2}         #y co-ordinates matrix of v grid
    xxc::Array{T,2}         #x co-ordinates matrix of c grid
    yyc::Array{T,2}         #y co-ordinates matrix of c grid
    σ::Vector{T}           #sigma levels
    ζ::Vector{T}            #reverse sigma levels
    quadrature_weights::Vector{T} #quadrature weights
end


#grid constructor 
function Grid(; 
                nx = 80,
                ny = 10,
                dx = 8000.0,
                dy = 8000.0,
                nσ = 4,
                x0 = 0.0,
                y0 = -40000.0,
                h_mask = trues(nx,ny),
                u_iszero = falses(nx+1,ny),
                v_iszero = falses(nx,ny+1))

    #check the sizes of inputs
    @assert size(h_mask)==(nx,ny);@assert h_mask == clip(h_mask)
    @assert size(u_iszero)==(nx+1,ny)
    @assert size(v_iszero)==(nx,ny+1)

    #map bit arrays to boolean
    h_mask = convert(Array{Bool,2}, h_mask)
    u_iszero = convert(Array{Bool,2}, u_iszero)
    v_iszero = convert(Array{Bool,2}, v_iszero)


    #compute grid co-ordinates
    xxh=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]; @assert size(xxh)==(nx,ny)
    yyh=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]; @assert size(yyh)==(nx,ny)

    xxu=[x0+(i-1.0)*dx for i=1:nx, j=1:ny]; @assert size(xxu)==(nx,ny)
    yyu=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]; @assert size(yyu)==(nx,ny)

    xxv=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]; @assert size(xxv)==(nx,ny)
    yyv=[y0+(j-1.0)*dy for i=1:nx, j=1:ny]; @assert size(yyv)==(nx,ny)

    xxc=[x0+i*dx for i=1:nx, j=1:ny]; @assert size(xxc)==(nx,ny)
    yyc=[y0+j*dy for i=1:nx, j=1:ny]; @assert size(yyc)==(nx,ny)
    
    #sigma grid info
    σ = collect(range(0.0,length=nσ,stop=1.0)); @assert length(σ) == nσ
    ζ = one(eltype(σ)) .- σ ; @assert length(ζ) == nσ
    quadrature_weights = [0.5;ones(nσ-2);0.5]/(nσ-1); @assert length(quadrature_weights) == nσ

    return Grid(nx,ny,nσ,dx,dy,x0,y0,h_mask,u_iszero,v_iszero,
                xxh,yyh,xxu,yyu,xxv,yyv,xxc,yyc,σ,ζ,quadrature_weights)
end



#Struct to hold model parameters.
#Format: fieldname::Type = default_value.
#T & N are type parameters, usually real numbers (e.g. Float64) and integers (e.g. Int64) respectively.

#bed_elevation::Array{T,2} = zeros(nx,ny); @assert size(bed_elevation)==(nx,ny)
#starting_thickness::Array{T,2} = zeros(nx,ny); @assert size(starting_thickness)==(nx,ny)

@with_kw struct Params{T <: Real, N <: Integer}
dt::T = 1.0
g::T = 9.81
density_ice::T = 918.0
density_ocean::T = 1028.0
gas_const=8.314;
sec_per_year::T = 3600*24*365.25
default_thickness::T = 100.
default_viscosity::T = 1.0e7
default_temperature::T = 265.700709
default_damage::T = 0.0
accumulation_rate::T = 0.0
basal_melt_rate::T = 0.0
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

#structure to store initial conditions
@with_kw struct InitialConditions{T <: Real}
    initial_thickness::Array{T,2} = zeros(10,10)
    initial_viscosity::Array{T,2} = zeros(10,10) #placeholder array
    initial_temperature::Array{T,2} = zeros(10,10)
    initial_damage::Array{T,2} = zeros(10,10)
end


#structure to hold the solver parameters
@with_kw struct SolverParams{T <: Real, N <: Integer}
    n_iter_viscosity::N = 2;  @assert n_iter_viscosity ==2
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

struct TimesteppingParams{T <: Real, N <: Integer, TO, C, P}
    n_iter0::N      #initial iteration number
    dt::T           #timestep
    end_time::T     #end time of this simulation
    t0::T           #start time of this simulation 
    chkpt_freq::T   #temporary checkpoint frequency
    pchkpt_freq::T  #permanent checkpoint frequency  
    n_iter_total::TO #total number of timesteps counting from zero
    n_iter_chkpt::C #number of iterations per temporary checkpoint
    n_iter_pchkpt::P#number of iterations per permanent checkpoint
end

function TimesteppingParams(;
                        n_iter0 = 0,
                        dt = 1.0,
                        end_time = 1.0,
                        t0 = nothing,
                        chkpt_freq = Inf,
                        pchkpt_freq = Inf)

    #(n_iter0 > 0) || ArgumentError("n_iter0 must be a positive number")

    #if n_iter0 > 0, check file exists and get start time, else throw error 


    #initialize t0 (really you should read start time from pickup file)
    t0 = n_iter0 > 0 ? n_iter0 * dt : 0 
    t0 = map(typeof(dt), t0)

    #compute number of timesteps (total and per checkpoint)
    end_time == Inf ? n_iter_total = Inf : n_iter_total  = round(Int, end_time/dt)
    chkpt_freq == Inf ? n_iter_chkpt = Inf : n_iter_chkpt  = round(Int, chkpt_freq/dt)
    pchkpt_freq == Inf ? n_iter_pchkpt = Inf : n_iter_pchkpt = round(Int, pchkpt_freq/dt)
    
    return TimesteppingParams(n_iter0, dt, end_time, t0, chkpt_freq, pchkpt_freq, n_iter_total, n_iter_chkpt, n_iter_pchkpt)
end


#mutable clock structure to store time info
mutable struct Clock{T <: Real, N <: Integer}
    n_iter::N
    time::T
end

#clock constructor
function Clock(;
                n_iter = 0,
                time = 0)
    return Clock(n_iter, time)
end

########### H grid explicit constructor ###########
#Struct to hold information on h-grid, located at cell centers.
#struct HGrid{T <: Real, N <: Integer}
#mask::Array{Bool,2}
#n::N
#crop::Diagonal{T,Array{T,1}}
#samp::SparseMatrixCSC{T,N} 
#spread::SparseMatrixCSC{T,N}
#b::Array{T,2} 
#h::Array{T,2} 
#s::Array{T,2} 
#dhdt::Array{T,2} 
#accumulation::Array{T,2}
#basal_melt::Array{T,2}
#haf::Array{T,2}
#grounded_fraction::Array{T,2}
#dsdh::Array{T,2}
#shelf_strain_rate::Array{T,2} 
#av_speed::Array{T,2}
#bed_speed::Array{T,2}
#weertman_c::Array{T,2}
#β::Array{T,2} 
#βeff::Array{T,2} 
#τbed::Array{T,2}
#ηav::Array{T,2}
#quad_f2::Array{T,2}
#dneghηav::Base.RefValue{Diagonal{T,Array{T,1}}} 
#dimplicit::Base.RefValue{Diagonal{T,Array{T,1}}}
#end

#H grid constructor
#function HGrid(grid::Grid{T,N}, params) where {T <: Real, N <: Integer}
#    #unpack grid size
#    nx = grid.nx
#    ny = grid.ny

#    #fill fields
#    mask = grid.h_mask
#    n = count(mask);                 @assert n == count(mask)
#    crop = Diagonal(float(mask[:])); @assert crop == Diagonal(float(mask[:]))
#    samp = crop[mask[:],:];          @assert samp == crop[mask[:],:]
#    spread = sparse(samp');          @assert spread == sparse(samp')
#    b = grid.bed_elevation;             @assert size(b)==(nx,ny)
#    h = params.starting_thickness;   @assert size(h)==(nx,ny)
#    s = zeros(nx,ny);                @assert size(s)==(nx,ny)
#    dhdt = zeros(nx,ny);             @assert size(dhdt)==(nx,ny)
#    accumulation = zeros(nx,ny);     @assert size(accumulation)==(nx,ny)
#    basal_melt = zeros(nx,ny);       @assert size(basal_melt)==(nx,ny)
#    haf = zeros(nx,ny);              @assert size(haf)==(nx,ny)
#    grounded_fraction = ones(nx,ny); @assert size(grounded_fraction)==(nx,ny)
#    dsdh = ones(nx,ny);              @assert size(dsdh)==(nx,ny)
#    shelf_strain_rate = zeros(nx,ny);@assert size(shelf_strain_rate)==(nx,ny)
#    av_speed = zeros(nx,ny);         @assert size(av_speed)==(nx,ny)
#    bed_speed = zeros(nx,ny);        @assert size(bed_speed)==(nx,ny)
#    weertman_c = zeros(nx,ny);       @assert size(weertman_c)==(nx,ny)
#    β = zeros(nx,ny);                @assert size(β)==(nx,ny)
#    βeff = zeros(nx,ny);             @assert size(βeff)==(nx,ny)
#    τbed = zeros(nx,ny);             @assert size(τbed)==(nx,ny)
#    ηav = fill(params.starting_viscosity,nx,ny)
#    quad_f2 = h./(3*ηav);           @assert size(quad_f2)==(nx,ny)
#    dneghηav = Ref(crop*Diagonal(zeros(nx*ny))*crop)
#    dimplicit = Ref(crop*Diagonal(zeros(nx*ny))*crop)
#
#    return HGrid(mask,n,crop,samp,spread,b,h,s,dhdt,accumulation,basal_melt,haf,
#        grounded_fraction,dsdh,shelf_strain_rate,av_speed,bed_speed,weertman_c,β,
#        βeff,τbed,ηav,quad_f2,dneghηav,dimplicit)
#end


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
grid::Grid{T,N}
params::Params{T,N}
timestepping_params::TimesteppingParams{T,N}
solver_params::SolverParams{T,N}
initial_conditions::InitialConditions{T}
gh::HGrid{T,N}
gu::UGrid{T,N}
gv::VGrid{T,N}
gc::CGrid{T,N}
g3d::SigmaGrid{T,N}   
wu::UWavelets{T,N}
wv::VWavelets{T,N}
clock::Clock{T,N}
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
function start(;
    grid = nothing, 
    bed_elevation = nothing,
    params = Params(),
    solver_params = SolverParams(),
    initial_conditions = InitialConditions(),
    timestepping_params = TimesteppingParams())

    #check that a grid and bed has been inputted
    ~(grid === nothing) || throw(ArgumentError("You must specify an input grid"))
    ~(bed_elevation === nothing) || throw(ArgumentError("You must input a bed elevation"))
    
    #check types
    #if a functional bed has been specified, convert to an array
    bed_array = zeros(grid.nx,grid.nx) #initialize a bed array
    try
    bed_array = get_bed_elevation(bed_elevation, grid)
    catch
    throw(ArgumentError("bed elevation must be of type function or array"))
    end
            #check the size of the bed
    #(Base.size(bed_array) = (grid.nx, grid.ny)) || throw(ArgumentError("Bed and grid sizes must be identical"))
    

    #Check initial conditions, and revert to default values if not
    initial_conditions = check_initial_conditions(initial_conditions, params, grid)

    #Define masks for points on h-, u-, v- and c-grids that lie in model domain.
    h_mask = grid.h_mask 
    u_mask = get_u_mask(h_mask)
    v_mask = get_v_mask(h_mask)
    c_mask = get_c_mask(h_mask)

    #Remove all points on u- and v-grids with homogenous Dirichlet conditions.
    u_mask[grid.u_iszero].=false
    v_mask[grid.v_iszero].=false

    #h-grid
    #gh=HGrid(grid, params) #uncomment if using the explicit constructor method
    gh=HGrid(
    x0=grid.x0,
    y0=grid.y0,
    nx=grid.nx,
    ny=grid.ny,
    dx=grid.dx,
    dy=grid.dy,
    mask=h_mask,
    b = bed_array,
    h = initial_conditions.initial_thickness,
    ηav = initial_conditions.initial_viscosity,
    )

    #u-grid
    gu=UGrid(
    x0=grid.x0,
    y0=grid.y0,
    nx=grid.nx+1,
    ny=grid.ny,
    dx=grid.dx,
    dy=grid.dy,
    mask=u_mask,
    levels=solver_params.levels
    )

    #v-grid
    gv=VGrid(
    x0=grid.x0,
    y0=grid.y0,
    nx=grid.nx,
    ny=grid.ny+1,
    dx=grid.dx,
    dy=grid.dy,
    mask=v_mask,
    levels=solver_params.levels
    )

    #c-grid
    gc=CGrid(
    x0=grid.x0,
    y0=grid.y0,
    nx=grid.nx-1,
    ny=grid.ny-1,
    dx=grid.dx,
    dy=grid.dy,
    mask=c_mask
    )

    #3D-grid
    g3d=SigmaGrid(
    nx=grid.nx,
    ny=grid.ny,
    nσ=grid.nσ,
    η = fill(params.default_viscosity,grid.nx,grid.ny,grid.nσ),
    θ = fill(params.default_temperature,grid.nx,grid.ny,grid.nσ),
    Φ = fill(params.default_damage,grid.nx,grid.ny,grid.nσ),
    glen_b = fill(glen_b(params.default_temperature,params.default_damage,params),grid.nx,grid.ny,grid.nσ)
    )

    #Wavelet-grid, u-component.
    wu=UWavelets(nx=grid.nx+1,ny=grid.ny,levels=params.levels)

    #Wavelet-grid, v-component.
    wv=VWavelets(nx=grid.nx,ny=grid.ny+1,levels=params.levels)

    #Default clock
    clock = Clock(n_iter = 0, time = 0.0)

    #Use type constructor to build initial state.
    wavi=State(grid,params,timestepping_params,solver_params,initial_conditions,gh,gu,gv,gc,g3d,wu,wv,clock)

    return wavi
end

"""
    check_initial_conditions(initial_conditions, params)

Check whether initial conditions have been specified. Default them to standard values if not
"""
function check_initial_conditions(initial_conditions, params, grid)
    if initial_conditions.initial_thickness == zeros(10,10)
        default_thickness = params.default_thickness
        println("Did not find a specified initial thickness, reverting to default value specified in params ($default_thickness m everywhere)")
        initial_conditions = @set initial_conditions.initial_thickness =  default_thickness*ones(grid.nx, grid.ny)
    end
    
    if initial_conditions.initial_viscosity == zeros(10,10)
        default_viscosity = params.default_viscosity
        println("Did not find a specified initial viscosity, reverting to default value specified in params ($default_viscosity Pa s everywhere)")
        initial_conditions = @set initial_conditions.initial_viscosity =  default_viscosity*ones(grid.nx, grid.ny)
    end

    if initial_conditions.initial_temperature == zeros(10,10)
        default_temperature = params.default_temperature
        println("Did not find a specified initial temperature, reverting to default value specified in params ($default_temperature K everywhere)")
        initial_conditions = @set initial_conditions.initial_temperature =  default_temperature*ones(grid.nx, grid.ny)
    end

    if initial_conditions.initial_damage == zeros(10,10)
        default_damage = params.default_damage
        println("Did not find a specified initial damage field, reverting to default value specified in params ($default_damage everywhere)")
        initial_conditions = @set initial_conditions.initial_damage =  default_damage*ones(grid.nx, grid.ny)
    end
    return initial_conditions
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
    update_clock!(wavi)
    return wavi
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
    simulation
    
"""
function simulation(; 
    grid = nothing,
    bed_elevation = nothing,
    params = nothing,
    solver_params = nothing,
    initial_conditions = nothing,
    timestepping_params = nothing)

    ~(timestepping_params === nothing) || error("Must pass timestepping params")

    if timestepping_params.n_iter0 == 0 #start a fresh run
        #check that grid and bed have been passed
        ~(grid === nothing) || error("Must pass a grid if starting a fresh run")
        ~(bed_elevation === nothing) || error("Must pass a bed (array or function) if starting a fresh run")

        println("Starting clean wavi simulation")


        #if no parameters have been passed, construct defaults
        if (params === nothing); params = Params(); end 
        if (solver_params === nothing); solver_params = SolverParams(); end 
        if (initial_conditions === nothing); initial_conditions = InitialConditions(); end #don't worry about strange defaults here, these will be picked up by start
        if (timestepping_params === nothing); timestepping_params = TimesteppingParams(); end 
        wavi = start(grid = grid, 
                    bed_elevation = bed_elevation,
                    params = params, 
                    solver_params = solver_params,
                    initial_conditions = initial_conditions,
                    timestepping_params = timestepping_params)

        #do the run
        chkpt_tag = "A" #initialize the checkpoint tag
        println("running simulation...")
        for i = 1:timestepping_params.n_iter_total
            run!(wavi)
            if mod(i,timestepping_params.n_iter_chkpt) == 0
                #output a temporary checkpoint
                fname = string("Chkpt",chkpt_tag, ".jld2")
                @save fname wavi
                chkpt_tag = (chkpt_tag == "A") ? "B" : "A"
                println("making temporary checkpoint at iteration number $(wavi.clock.n_iter)")
            end
            if mod(i,timestepping_params.n_iter_pchkpt) == 0
                #output a permanent checkpoint
                n_iter_string =  lpad(wavi.clock.n_iter, 10, "0"); #filename as a string with 10 digits
                fname = string("PChkpt_",n_iter_string, ".jld2")
                @save fname wavi
                println("making permanent checkpoint at iteration number $(wavi.clock.n_iter)")

            end
        end
        

    else #look for a pickup
        n_iter_string =  lpad(timestepping_params.n_iter0, 10, "0"); #filename as a string with 10 digits
        try 
            @load string("PChkpt_",n_iter_string, ".jld2") wavi
            println("Pickup successful")
        catch 
            println("Pickup error, terminating run")
        end

        #update the parameters of those that have been specified the flag is specified

        #continue with the run
        chkpt_tag = "A" #initialize the checkpoint tag
        for i = (timestepping_params.n_iter0+1):timestepping_params.n_iter_total
            run!(wavi)
            if mod(i,timestepping_params.n_iter_chkpt) == 0
                #output a temporary checkpoint
                println("making temporary checkpoint at iteration number $(wavi.clock.n_iter)")
                fname = string("Chkpt",chkpt_tag, ".jld2")
                @save fname wavi
                chkpt_tag = (chkpt_tag == "A") ? "B" : "A"
            end
            if mod(i,timestepping_params.n_iter_pchkpt) ==0 
                #output a permanent checkpoint
                println("making permanent checkpoint at iteration number $(wavi.clock.n_iter)")
                n_iter_string =  lpad(wavi.clock.n_iter, 10, "0"); #filename as a string with 10 digits
                fname = string("PChkpt_",n_iter_string, ".jld2")
                @save fname wavi

            end
        end
    end
    return wavi
end

# Utility functions
function get_bed_elevation(bed_elevation::F, grid) where (F <: Function)
    bed_array = bed_elevation.(grid.xxh, grid.yyh)
    return bed_array
end

function get_bed_elevation(bed_elevation::Array{T,2}, grid) where (T <: Real)
    bed_array = bed_elevation
    return bed_array
end


"""
     get_glx(wavi)

Return the grounding line in the form x = x(y). Assumes each y-row has at least
one grid point where ice grounded and one where ice floating.
"""
function get_glx(wavi)
      glmask=diff(sign.(wavi.gh.haf),dims=1).==-2 #calculate where sign of height above floating passes thru zero
      glx1=wavi.grid.xxh[1:end-1,:][glmask] #x co-ordiates upstream of grounding line
      glx2=wavi.grid.xxh[2:end,:][glmask] #x co-ordinates immediately downstream
      haf1=wavi.gh.haf[1:end-1,:][glmask] #Height above floatation immediately upstream of gl
      haf2=wavi.gh.haf[2:end,:][glmask]
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
    update_clock!(wavi::AbstractModel)

Update the clock

"""
function update_clock!(wavi::AbstractModel)
    @unpack clock, timestepping_params = wavi
    clock.n_iter += 1
    clock.time += timestepping_params.dt
    return wavi
end


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
    @unpack gu,gv,wu,wv,params,solver_params=wavi

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
        converged = rel_resid < solver_params.tol_picard

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
    @unpack gu,gv,params,solver_params=wavi
    m,n=size(op)
    @assert m == n == gu.n + gv.n
    op_diag=zeros(eltype(op),n)
    sm=solver_params.stencil_margin
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

    @unpack gu,gv,wu,wv,params,solver_params=wavi

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
            op_diag=op_diag, nsmooth=solver_params.nsmooth, tol_coarse = solver_params.tol_coarse,
            maxiter_coarse = solver_params.maxiter_coarse, smoother_omega=solver_params.smoother_omega)

    return p
end
"""
    LinearAlgebra.ldiv!(x::AbstractVecOrMat{T}, p::AbstractPreconditioner{T,N}, b::AbstractVecOrMat{T}) where {T,N}

Overload LinearAlgebra.ldiv! function so that the bespoke preconditioner is deployed in calls to the
conjugate gradient method if p has type  <: AbstractPreconditioner.

"""
function ldiv!(x::AbstractVecOrMat{T}, p::AbstractPreconditioner{T,N}, b::AbstractVecOrMat{T}) where {T,N}
    precondition!(x, p, b)
end
"""
    precondition!(x, p, b)

Apply wavelet-based multigrid preconditioner using information stored in p.
"""
function precondition!(x, p, b)
    @unpack op,op_diag,nsmooth,sweep,sweep_order,smoother_omega,restrict,
            prolong,op_coarse,correction_coarse,tol_coarse,maxiter_coarse = p


    n=size(op,1)

    # Multigrid smooth
    x .= gauss_seidel_smoother!(x, op, b; iters = nsmooth, op_diag=op_diag,
                                sweep=sweep, sweep_order=sweep_order, smoother_omega = smoother_omega)

    resid=b-op*x;

    # Multigrid restriction
    b_coarse=restrict*resid

    # Multigrid solve for correction
    cg!(correction_coarse, op_coarse, b_coarse; reltol = tol_coarse, maxiter = maxiter_coarse)

    # Multigrid prolongation
    x .= x .+ prolong*correction_coarse

    # Multigrid smooth
    x .= gauss_seidel_smoother!(x, op, b; iters = nsmooth, op_diag=op_diag,
                                sweep=sweep, sweep_order=reverse(sweep_order), smoother_omega = smoother_omega)

    return x
end
"""
    gauss_seidel_smoother!(x, op, b;
                                iters=5,
                                op_diag=diag(sparse(op)),
                                sweep=(1:size(op,1)),
                                sweep_order=unique(sweep),
                                smoother_omega=1.0)
Apply smoother used in multigrid preconditioner.
"""
function gauss_seidel_smoother!(x, op, b;
                                iters=5,
                                op_diag=diag(sparse(op)),
                                sweep=(1:size(op,1)),
                                sweep_order=unique(sweep),
                                smoother_omega=1.0)
    resid=b-op*x
    for i = 1:iters
        for j = sweep_order
              idx = sweep .== j
              x[idx] .= x[idx] .+ smoother_omega .* resid[idx]./op_diag[idx]
              resid .= b .- op*x
        end
    end
    return x
end
"""
    update_surface_elevation!(wavi::AbstractModel)

Adjust surface elevation to hydrostatic equilibrium.
"""
function update_surface_elevation!(wavi::AbstractModel)
    @unpack gh,params=wavi
    gh.s[gh.mask] .= max.(gh.b[gh.mask]+gh.h[gh.mask],
                          params.sea_level_wrt_geoid .+ gh.h[gh.mask]*(1-params.density_ice./params.density_ocean))
    return wavi
end
"""
    update_dhdt!(wavi::AbstractModel)

Evaluate rate of change of thickness using mass conservation.
"""
function update_dhdt!(wavi::AbstractModel)
    @unpack gh,gu,gv=wavi
    gh.dhdt[gh.mask].=gh.samp*(gh.accumulation[:] .- gh.basal_melt[:] .-
             (  (gu.∂x*(gu.crop*(gu.h[:].*gu.u[:]))) .+ (gv.∂y*(gv.crop*(gv.h[:].*gv.v[:]))) ) )
    return wavi
end
"""
    update_thickness!(wavi::AbstractModel)

Update thickness using rate of change of thickness and apply minimum thickness constraint.
"""
function update_thickness!(wavi::AbstractModel)
    @unpack gh,gu,gv,params=wavi
    onesvec=ones(gh.nx*gh.ny)
    gh.h[gh.mask] .= gh.h[gh.mask] .+ max.(params.minimum_thickness .- gh.h[gh.mask],params.dt*gh.dhdt[gh.mask])
    return wavi
end
"""
    update_height_above_floatation!(wavi::AbstractModel)

Update height above floatation. Zero value is used to define location of grounding line.
"""
function update_height_above_floatation!(wavi::AbstractModel)
    @unpack gh,params=wavi
    gh.haf .= height_above_floatation.(gh.h,gh.b,Ref(params))
    return wavi
end
"""
    update_grounded_fraction_on_huv_grids!(wavi::AbstractModel)

Update grounded area fraction on h-, u-, and v-grids for use in subgrid parameterisation.
"""
function update_grounded_fraction_on_huv_grids!(wavi::AbstractModel)
    @unpack gh,gu,gv = wavi
    (gfh,gfu,gfv)=pos_fraction(gh.haf;mask=gh.mask)
    gh.grounded_fraction[:] .= gfh[:]
    gu.grounded_fraction[:] .= gfu[:]
    gv.grounded_fraction[:] .= gfv[:]
    return wavi
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
"""
    update_dsdh!(wavi::AbstractModel)

Compute change of surface elevation per unit thickness change, accounting for hydrostatic adjustment.
"""
function update_dsdh!(wavi::AbstractModel)
    @unpack gh,gu,gv,params = wavi
    gh.dsdh .= (1.0 - params.density_ice./params.density_ocean) .+
           (params.density_ice./params.density_ocean).*gh.grounded_fraction;
    return wavi
end
"""
    update_weertman_c!(wavi::AbstractModel)

Update coefficient used in the sliding law to account for migration of grounding line.
"""
function update_weertman_c!(wavi::AbstractModel)
    @unpack gh,params = wavi
    gh.weertman_c .= params.weertman_c .* gh.grounded_fraction
    return wavi
end
"""
    update_accumulation_rate!(wavi::AbstractModel)

Update the accumulation rate.
"""
function update_accumulation_rate!(wavi::AbstractModel)
    @unpack gh,params = wavi
    gh.accumulation .= params.accumulation_rate
    return wavi
end
"""
    update_basal_melt!(wavi::AbstractModel)

Update the basal melt rate.
"""
function update_basal_melt!(wavi::AbstractModel)
    @unpack gh, params = wavi
    gh.basal_melt .= params.basal_melt_rate
    return wavi
end
"""
    update_βeff!(wavi::AbstractModel)

Compute the effective drag coefficient.
"""
function update_βeff!(wavi::AbstractModel)
    @unpack gh=wavi
    gh.βeff .= gh.β ./ (1.0 .+ gh.quad_f2 .* gh.β)
    return wavi
end
"""
    update_βeff_on_uv_grids!(wavi::AbstractModel)

Interpolate the effective drag coefficient onto u- and v-grids, accounting for grounded fraction.
"""
function update_βeff_on_uv_grids!(wavi::AbstractModel)
    @unpack gh,gu,gv=wavi
    @assert eltype(gh.grounded_fraction)==eltype(gh.βeff)

    T=eltype(gh.grounded_fraction)

    onesvec=ones(T,gh.nx*gh.ny)
    gu.βeff[gu.mask].=(gu.samp*(gu.cent'*(gh.crop*gh.βeff[:])))./(gu.samp*(gu.cent'*(gh.crop*onesvec)))
    ipolgfu=zeros(T,gu.nx,gu.ny);
    ipolgfu[gu.mask].=(gu.samp*(gu.cent'*(gh.crop*gh.grounded_fraction[:])))./(gu.samp*(gu.cent'*(gh.crop*onesvec)))
    gu.βeff[ipolgfu .> zero(T)] .= gu.βeff[ipolgfu .> zero(T)].*gu.grounded_fraction[ipolgfu .> zero(T)]./
                                                        ipolgfu[ipolgfu .> zero(T)]

    gv.βeff[gv.mask].=(gv.samp*(gv.cent'*(gh.crop*gh.βeff[:])))./(gv.samp*(gv.cent'*(gh.crop*onesvec)))
    ipolgfv=zeros(T,gv.nx,gv.ny);
    ipolgfv[gv.mask].=(gv.samp*(gv.cent'*(gh.crop*gh.grounded_fraction[:])))./(gv.samp*(gv.cent'*(gh.crop*onesvec)))
    gv.βeff[ipolgfv .> zero(T)] .= gv.βeff[ipolgfv .> zero(T)].*gv.grounded_fraction[ipolgfv .> zero(T)]./
                                                 ipolgfv[ipolgfv .> zero(T)];

    return wavi
end
"""
    set_velocities!(wavi::AbstractModel,x)

Set velocities to particular values. Input vector x represents stacked u and v components at valid grid points.
"""
function set_velocities!(wavi::AbstractModel,x)
    @unpack gh,gu,gv,gc=wavi
    gu.u[:] .= gu.spread*x[1:gu.n]
    gv.v[:] .= gv.spread*x[(gu.n+1):(gu.n+gv.n)]
    return wavi
end
"""
    update_shelf_strain_rate!(wavi::AbstractModel)

Find the effective strain rate for 'ice shelf' parts of strain rate tensor, neglecting all vertical shear.
"""
function update_shelf_strain_rate!(wavi::AbstractModel)
    @unpack gh,gu,gv,gc=wavi
    gh.shelf_strain_rate[:] .= sqrt.( (gh.crop*(gu.∂x*(gu.crop*gu.u[:]))).^2 .+
                                      (gh.crop*(gv.∂y*(gv.crop*gv.v[:]))).^2 .+
                                (gh.crop*(gu.∂x*(gu.crop*gu.u[:]))).*(gh.crop*(gv.∂y*(gv.crop*gv.v[:]))) .+
                       0.25*(gh.crop*(gc.cent*(gc.crop*( gu.∂y*(gu.crop*gu.u[:]) .+ gv.∂x*(gv.crop*gv.v[:]) )))).^2  )
    return wavi
end
"""
    update_av_speed!(wavi::AbstractModel)

Find the depth-averaged speed on the h-grid using components on u- and v- grids
"""
function update_av_speed!(wavi::AbstractModel)
    @unpack gh,gu,gv=wavi
    gh.av_speed[:] .= sqrt.( (gh.crop*(gu.cent*(gu.crop*gu.u[:]))).^2 .+ (gh.crop*(gv.cent*(gv.crop*gv.v[:]))).^2 )
    return wavi
end
"""
    update_bed_speed!(wavi::AbstractModel)

Find the sliding speed at the bed on the h-grid using the average speed.
"""
function update_bed_speed!(wavi::AbstractModel)
    @unpack gh=wavi
    gh.bed_speed .= gh.av_speed ./ (1.0 .+ gh.quad_f2 .* gh.β)
    return wavi
end
"""
    update_β!(wavi::AbstractModel)

Find the drag coefficient at the bed using the sliding law.
"""
function update_β!(wavi::AbstractModel)
    @unpack gh,params=wavi
    gh.β .= gh.weertman_c .* ( sqrt.(gh.bed_speed.^2 .+  params.weertman_reg_speed^2 ) ).^(1.0/params.weertman_m - 1.0)
    return wavi
end
"""
    update_basal_drag!(wavi::AbstractModel)

Find the shear stress at the bed.
"""
function update_basal_drag!(wavi::AbstractModel)
    @unpack gh=wavi
    gh.τbed .= gh.β .* gh.bed_speed
    return wavi
end
"""
    inner_update_viscosity!(wavi::AbstractModel)

Inner update to iteratively refine viscosity on the 3d grid at all sigma levels.
"""
function inner_update_viscosity!(wavi::AbstractModel)
    @unpack gh,g3d,params,solver_params=wavi
    for k=1:g3d.nσ
        for j=1:g3d.ny
            for i=1:g3d.nx
                if gh.mask[i,j]
                    for iter=1:solver_params.n_iter_viscosity
                        g3d.η[i,j,k] = 0.5 * g3d.glen_b[i,j,k] * (
                                                   sqrt(    gh.shelf_strain_rate[i,j]^2 +
                                                            0.25*(gh.τbed[i,j]*g3d.ζ[k]/g3d.η[i,j,k])^2 +
                                                            params.glen_reg_strain_rate^2   )
                                                                 )^(1.0/params.glen_n - 1.0)
                    end
                end
            end
        end
    end
    return wavi
end
"""
    update_av_viscosity!(wavi::AbstractModel)

Use quadrature to compute the depth averaged viscosity.
"""
function update_av_viscosity!(wavi::AbstractModel)
    @unpack gh,g3d=wavi
    gh.ηav .= zero(gh.ηav)
    for k=1:g3d.nσ
       for j = 1:g3d.ny
          for i = 1:g3d.nx
                gh.ηav[i,j] += g3d.quadrature_weights[k] * g3d.η[i,j,k]
          end
       end
    end
    return wavi
end
"""
    update_av_viscosity!(wavi::AbstractModel)

Use quadrature to compute f2 function, used to relate average velocities to basal velocities.
"""
function update_quadrature_f2!(wavi::AbstractModel)
    @unpack gh,g3d=wavi
    gh.quad_f2 .= zero(gh.quad_f2)
    for k=1:g3d.nσ
       for j = 1:g3d.ny
          for i = 1:g3d.nx
                gh.quad_f2[i,j] += g3d.quadrature_weights[k]*gh.h[i,j]*(g3d.ζ[k])^2/g3d.η[i,j,k]
          end
       end
    end
    return wavi
end
"""
    update_rheological_operators!(wavi::AbstractModel)

Precompute various diagonal matrices used in defining the momentum operator.
"""
function update_rheological_operators!(wavi::AbstractModel)
    @unpack gh,gu,gv,params=wavi
    gh.dneghηav[] .= gh.crop*Diagonal(-gh.h[:].*gh.ηav[:])*gh.crop
    gu.dnegβeff[] .= gu.crop*Diagonal(-gu.βeff[:])*gu.crop
    gv.dnegβeff[] .= gv.crop*Diagonal(-gv.βeff[:])*gv.crop
    gh.dimplicit[] .= gh.crop*Diagonal(-params.density_ice * params.g * params.dt * gh.dsdh[:])*gh.crop
    return wavi
end

"""
    opvec(wavi::AbstractModel,vec::AbstractVector)

Function to multiply a vector by the momentum operator.
"""
function opvec(wavi::AbstractModel,vec::AbstractVector)
    @unpack gh,gu,gv,gc=wavi
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
    restrictvec(wavi::AbstractModel,vec::AbstractVector)

Function to restrict a vector from the fine grid to the coarse grid, used in multigrid preconditioner.
"""
function restrictvec(wavi::AbstractModel,vec::AbstractVector)
    @unpack wu,wv,gu,gv=wavi
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
    prolongvec(wavi::AbstractModel,waveletvec::AbstractVector)

Function to prolong a vector from the coarse grid to the fine grid, used in multigrid preconditioner.
"""
function prolongvec(wavi::AbstractModel,waveletvec::AbstractVector)
    @unpack wu,wv,gu,gv=wavi
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
    update_wavelets!(wavi::AbstractModel)

Compute wavelet transform of velocities to define the coarse grid used in multigrid preconditioner.
"""
function update_wavelets!(wavi::AbstractModel)
    @unpack wu,wv,gu,gv,params,solver_params=wavi

    wu.wavelets[:] .= gu.dwt*(gu.crop*gu.u[:])
    wv.wavelets[:] .= gv.dwt*(gv.crop*gv.v[:])

    wu.mask .= (abs.(wu.wavelets) .>= solver_params.wavelet_threshold)
    wv.mask .= (abs.(wv.wavelets) .>= solver_params.wavelet_threshold)

    wu.n[] = count(wu.mask)
    wv.n[] = count(wv.mask)

    wu.crop[] .= Diagonal(float(wu.mask[:]))
    wv.crop[] .= Diagonal(float(wv.mask[:]))

    wu.samp[] = wu.crop[][wu.mask[:],:]
    wv.samp[] = wv.crop[][wv.mask[:],:]

    wu.spread[] = sparse(wu.samp[]')
    wv.spread[] = sparse(wv.samp[]')

    return wavi
end
"""
    wavelet_matrix(n,levels,direction)

Compute matrix used to apply wavelet transform in one dimension.
"""
function wavelet_matrix(n,levels,direction);

    identity=spdiagm(n,n, 0 => ones(n))

    isforward = lowercase(direction) == "forward"
    isreverse = lowercase(direction) == "reverse"

    if isforward
        fwav=identity
    elseif isreverse
        rwav=identity
    else
        error("Direction of wavelet transform not defined properly")
    end

    step=1

    for j=levels:-1:0

        step=step*2

        ixodd=1:step:(n-div(step,2))
        ixeven=ixodd.+div(step,2)
        ixothers=setdiff(1:n, union(ixodd,ixeven))

        odds=identity[ixodd,:]
        evens=identity[ixeven,:]
        others=identity[ixothers,:]

        if isforward
            #Forward Haar wavelet lifting scheme
            fwav = (others'*others+0.5*(evens'*(evens - odds) + odds'*(odds + evens)))*fwav
        elseif isreverse
            #Reverse Haar wavelet lifting scheme
            rwav = rwav*(others'*others + evens'*(evens + odds) + odds'*(odds - evens))
        else
            error("Direction of wavelet transform not defined properly")
        end

    end

    if isforward
            #Forward Haar wavelet lifting scheme
            return fwav
    elseif isreverse
            #Reverse Haar wavelet lifting scheme
            return rwav
    else
            error("Direction of wavelet transform not defined properly")
    end

end
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
