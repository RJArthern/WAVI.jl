export ThreadedSpec

using LinearAlgebra
using Parameters

using WAVI
using WAVI.Processes
import WAVI: AbstractGrid, AbstractModel, AbstractSpec
import WAVI.Processes: update_preconditioner!, precondition!

include("SchwarzDecomposition/SchwarzDecomposition.jl")

"""
Struct to represent the shared memory parallel specification of a model.

"""
@with_kw struct ThreadedSpec{T<:Real, N<:Integer} <: AbstractSpec
    ngridsx::N = 1 
    ngridsy::N = 1
    overlap::N = 1
    niterations::N = 0
    damping::T = 0.0
    schwarzModelArray::Array{AbstractModel,2} = Array{AbstractModel,2}(undef,ngridsx,ngridsy)
end

function update_preconditioner!(model::AbstractModel, spec::ThreadedSpec)
    @unpack ngridsx, ngridsy, overlap = spec
    n_patch = ngridsx * ngridsy
    if !isassigned(spec.schwarzModelArray, 1)
        @info "ThreadedSpec: building $(n_patch) Schwarz subdomain models"
    end
    @debug "Spawning $(n_patch) threads for preconditioning"

    @sync for igrid = 1:ngridsx
        for jgrid = 1:ngridsy
            Threads.@spawn begin
                spec.schwarzModelArray[igrid,jgrid] = schwarzModel(model;
                                                                         igrid=igrid,
                                                                         jgrid=jgrid,
                                                                         ngridsx=ngridsx,
                                                                         ngridsy=ngridsy, 
                                                                         overlap=overlap)
            end
        end
    end
    return model
end

update_preconditioner!(model::AbstractModel{T,N,S}) where {T,N,S<:ThreadedSpec} =
    update_preconditioner!(model, model.spec)


"""
precondition!(model::AbstractModel, spec::ThreadedSpec)

Apply an Additive Schwarz preconditioner with Partition-of-Unity (AS-PoU) weights
using shared memory parallelism.

# Method
For each Schwarz iteration:
1. **Restrict**: copy the current global velocities into each subdomain (overlapping regions).
2. **Local solve**: call `update_state!` on each subdomain independently (in parallel threads).
3. **Prolong (AS-PoU)**: accumulate weighted subdomain solutions back to the global domain.
"""
function precondition!(model::AbstractModel, spec::ThreadedSpec)
    @unpack ngridsx, ngridsy, overlap, niterations, schwarzModelArray, damping = spec
    @unpack solver_params = model

    @debug "Preconditioning across the $(ngridsx * ngridsy) threads"
    x = get_start_guess(model)  
    op = get_op(model)
    b = get_rhs(model)
    resid = get_resid(x,op,b)
    set_residual!(model,resid)
    rel_resid = norm(resid)/norm(b)
    converged = rel_resid < solver_params.tol_picard

    if ! converged
        for iteration = 1:niterations
            @debug "Schwarz iteration $iteration"
            @sync for igrid = 1:ngridsx
                for jgrid = 1:ngridsy
                    Threads.@spawn begin                
                        model_g = schwarzModelArray[igrid,jgrid]

                        schwarzRestrictVelocities!(
                            model_g::AbstractModel,
                            model::AbstractModel;
                            igrid=igrid,
                            jgrid=jgrid,
                            ngridsx=ngridsx,
                            ngridsy=ngridsy,
                            overlap=overlap)
                    end
                end
            end
            
            @sync for igrid = 1:ngridsx
                for jgrid = 1:ngridsy
                    Threads.@spawn begin
                        model_g = schwarzModelArray[igrid,jgrid]
                        update_state!(model_g)
                    end
                end
            end

            model.fields.gu.u[:,:] .= damping .* model.fields.gu.u
            model.fields.gv.v[:,:] .= damping .* model.fields.gv.v
            
            threadLock=ReentrantLock()
            @sync for igrid = 1:ngridsx
                for jgrid = 1:ngridsy
                    Threads.@spawn begin
                        model_g = schwarzModelArray[igrid,jgrid]

                        lock(threadLock)
                        try
                            schwarzProlongVelocities!(
                                model::AbstractModel,
                                model_g::AbstractModel;
                                igrid=igrid,
                                jgrid=jgrid,
                                ngridsx=ngridsx,
                                ngridsy=ngridsy,
                                overlap=overlap,
                                damping=damping)
                        finally
                            unlock(threadLock)
                        end
                    end                    
                end
            end

        end

    end
    return converged, rel_resid
end

precondition!(model::AbstractModel{<:Any, <:Any, <:ThreadedSpec}) =
    precondition!(model, model.spec)
