export SheetOnlyGlaDS

struct SheetOnlyGlaDS{T <: Real, N  <: Integer} <: AbstractBasalHydrology
    # accelerated pseudo-transient (APT) iterative loop - GlaDS numerics
    ncheck      :: N
    maxiter     :: N
    etol        :: T
    do_autotune :: Bool
    n_autotune  :: N
    ɛ           :: T

    # GlaDS physics
    k     :: T
    α     :: T
    β     :: T
    ϕ_fixed :: T
    Ã     :: T
    
    # water thickness - cavity opening (w)
    h_r :: T
    l_r :: T

    min_effective_pressure :: T
end

"""
SheetOnlyGlaDS(; <kwargs>)


Keyword arguments
=================
- ncheck                    : 
- maxiter                   : 
- etol                      : 
- do_autotune               : 
- n_autotune                :
- ɛ                         : avoids division by zero
- k                         : sheet conductivity
- α                         : 1st turbulent flow exponent
- β                         : 2nd turbulent flow exponent
- ϕ_fixed                   : hydraulic potential constant for fixed grid cells
- Ã                         : rheological constant of ice multiplied by an order-one geometrical factor
- h_r                       : bedrock bump height     [m]
- l_r                       : bedrock bump wavelength [m]
- min_effective_pressure    : minimum effective pressure (Pa)
"""

function SheetOnlyGlaDS(;
                    # APT iterative loop - GlaDS numerics
                    ncheck      = 10*128,
                    maxiter     = 10*128^2,
                    etol        = 1e-8,
                    do_autotune = true,
                    n_autotune  = 10,
                    ɛ           = 1e-9,

                    # GlaDS physics
                    k     = 0.005,    # Sheet conductivity
                    α     = 5 / 4,    # 1st turbulent flow exponent (GlaDS)
                    β     = 3 / 2,    # 2nd turbulent flow exponent
                    ϕ_fixed = 0.,     # hydraulic potential for fixed grid cells
                    Ã     = 2.5e-25,  # rheological constant of ice multiplied by an order-one geometrical factor
                    
                    # water thickness - cavity opening (w)
                    h_r = 0.1,        # bedrock bump height     [m]
                    l_r = 2.0,        # bedrock bump wavelength [m]

                    min_effective_pressure = 1.0e4
                    )

    return SheetOnlyGlaDS(
                    ncheck,
                    maxiter,
                    etol,
                    do_autotune,
                    n_autotune,
                    ɛ,
                    k,
                    α,
                    β,
                    ϕ_fixed,
                    Ã,
                    h_r,
                    l_r,
                    min_effective_pressure)
end

@views function residual!(r, ϕ, ∂ₓϕ, ∂ᵧϕ, qx, qy, D, ρigH, ρwgB, Neff, h, k, α, β, ϕ_fixed, u_b, h_r, l_r, Ã, glens_n, dx, dy, m, floatMask, domainMask)
    # p-Laplacian
    # q = -kh^α|∇ϕ|^(β-2)∇ϕ
    # where D := kh^α|∇ϕ|^(β-2)
    # boundary conditions
    @. ϕ *= domainMask.gh * (!floatMask)
    @. ϕ[domainMask.gphi] = ϕ_fixed # only required for GlaDS verification test

    # effective pressure
    @. Neff = ρigH + ρwgB - ϕ
    @. Neff *= (!floatMask)

    # compute hydraulic potential gradients
    # part c
    @. ∂ₓϕ.gu = (ϕ[2:end, :] - ϕ[1:end-1, :]) / dx
    @. ∂ᵧϕ.gv = (ϕ[:, 2:end] - ϕ[:, 1:end-1]) / dy
    
    # part uv
    @. ∂ₓϕ.gv[2:end-1, :] = (ϕ[3:end, 1:end-1] - ϕ[1:end-2, 1:end-1] + ϕ[3:end, 2:end] - ϕ[1:end-2, 2:end]) / (4*dx)
    @. ∂ᵧϕ.gu[:, 2:end-1] = (ϕ[1:end-1, 3:end] - ϕ[1:end-1, 1:end-2] + ϕ[2:end, 3:end] - ϕ[2:end, 1:end-2]) / (4*dy)
    
    # compute diffusivity
    ∇ϕr = 1e-12
    @. D.gu = (k * (h.gu)^α) * abs(∂ₓϕ.gu^2 + ∂ᵧϕ.gu^2 + ∇ϕr)^((β - 2) / 2)
    @. D.gv = (k * (h.gv)^α) * abs(∂ₓϕ.gv^2 + ∂ᵧϕ.gv^2 + ∇ϕr)^((β - 2) / 2)

    # update upper and lower part of each cell
    @. qx[2:end-1, :] = -D.gu * ∂ₓϕ.gu
    @. qx *= domainMask.gu

    # update left and right part of each cell
    @. qy[:, 2:end-1] = -D.gv * ∂ᵧϕ.gv
    @. qy *= domainMask.gv

    # compute residual
    # r = m - ∇·q + v - w
    @. r = -((qx[2:end, :] - qx[1:end-1, :]) / dx + (qy[:, 2:end] - qy[:, 1:end-1]) / dy) + m +
             (Ã * h.gh * abs(Neff)^(glens_n - 1) * (Neff)) -
             (max(u_b / l_r * (h_r - h.gh), 0.0))
    @. r *= domainMask.gh * (!floatMask)
    @. r[domainMask.gphi] *= 0. # only required for GlaDS verification test

    return
end

@views function jvp!(v, r̄, r, ϕ, ϕ̄, ∂ₓϕ, ∂ₓϕ̄, ∂ᵧϕ, ∂ᵧϕ̄, qx, qx̄, qy, qȳ, D, D̄, ρigH, ρigH̄, ρwgB, ρwgB̄, Neff, N̄eff, h, h̄, k, α, β, ϕ_fixed, u_b, u_b̄, h_r, l_r, Ã, glens_n, dx, dy, m, m̄, floatMask, domainMask)
    # Set the input vector for cell and vertex components
    @. ϕ̄  = v

    # Use Enzyme's forward mode AD to compute the JVP
    Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Forward),
                    residual!,
                    DuplicatedNoNeed(r, r̄),
                    DuplicatedNoNeed(ϕ, ϕ̄),
                    DuplicatedNoNeed(∂ₓϕ, ∂ₓϕ̄),
                    DuplicatedNoNeed(∂ᵧϕ, ∂ᵧϕ̄),
                    DuplicatedNoNeed(qx, qx̄),
                    DuplicatedNoNeed(qy, qȳ),
                    DuplicatedNoNeed(D, D̄),
                    DuplicatedNoNeed(ρigH, ρigH̄),
                    DuplicatedNoNeed(ρwgB, ρwgB̄),
                    DuplicatedNoNeed(Neff, N̄eff),
                    DuplicatedNoNeed(h, h̄),
                    Const(k),
                    Const(α),
                    Const(β),
                    Const(ϕ_fixed),
                    Duplicated(u_b,u_b̄),
                    Const(h_r),
                    Const(l_r),
                    Const(Ã),
                    Const(glens_n),
                    Const(dx),
                    Const(dy),
                    DuplicatedNoNeed(m, m̄),
                    Const(floatMask),
                    Const(domainMask))
    return nothing
end

"""
Extract diag(Jϕ) from J block.
- The coloring is based on diagonal block J of the Jacobian

Returns one tuple of indices containing two colors:
- I: Two groups of center grid indices for diag(J)
"""
function coloring(nx, ny)
    Nx_selected = LinearIndices((nx, ny))

    # indices for coloring of center grid
    m = iseven.((1:nx) .+ (1:ny)')
    I = (vec(Nx_selected[m]), vec(Nx_selected[.!m]))

    return I
end

"""
            update_basal_water_thickness_effective_pressure!(basal_hydrology::SheetOnlyGlaDS, model::AbstractModel; update_basal_water_thickness::Bool = true)
    
use the sheets only version of GlaDS to calculate the basal water thickness, hydraulic potential at the bed, and effective pressure
"""

function update_basal_water_thickness_effective_pressure!(basal_hydrology::SheetOnlyGlaDS{T}, model::AbstractModel; update_basal_water_thickness::Bool = true) where {T}
    # derived parameters
    @unpack ncheck, maxiter, etol, do_autotune, n_autotune, ɛ, k, α, β, ϕ_fixed, Ã, h_r, l_r = basal_hydrology
    @unpack nx, ny, dx, dy = model.grid
    @unpack gh, gu, gv = model.fields
    @unpack params = model
    
    # arrays
    ∂ₓϕ  = (gu=zeros(T,nx - 1, ny), gv=zeros(T,nx, ny - 1))
    ∂ᵧϕ  = (gu=zeros(T,nx - 1, ny), gv=zeros(T,nx, ny -1))
    qx   = zeros(T,nx + 1, ny)
    qy   = zeros(T,nx, ny + 1)
    r    = zeros(T,nx, ny)
    z    = zeros(T,nx, ny)
    z0   = zeros(T,nx, ny)
    s    = zeros(T,nx, ny)
    P    = zeros(T,nx, ny)
    D    = (gu=zeros(T,nx - 1, ny), gv=zeros(T,nx, ny - 1))
    ρigH = zeros(T,nx, ny)
    ρwgB = zeros(T,nx, ny) # elevation potential ϕ_m
    Neff = zeros(T,nx, ny) # effective pressure
    floatMask = falses(nx, ny)
    domainMask = (gh=gh.mask, gu=gu.mask, gv=gv.mask, gphi = gh.hyd_potential_isfixed)

    m    = gh.basal_melt./params.sec_per_year # m/s
    u_b  = gh.bed_speed./params.sec_per_year # m/s
    ϕ    = gh.hydraulic_potential_b

    # calculate potentials
    @. ρigH = params.density_ice * params.g * gh.h
    @. ρwgB = params.density_freshwater * params.g * gh.b

    # calculate floating area
    floatMask .= gh.grounded_fraction .< 1.0

    # set water sheet thickness
    h    = (gh=gh.basal_water_thickness, gu=zeros(T,nx - 1, ny), gv=zeros(T,nx, ny - 1)) # Water sheet thickness (shall compute initial thickness)
    @. h.gh *= domainMask.gh * (!floatMask)
    @. h.gu = abs(h.gh[2:end, :] + h.gh[1:end-1, :]) * 0.5
    @. h.gv = abs(h.gh[:, 2:end] + h.gh[:, 1:end-1]) * 0.5
    @. h.gu *= domainMask.gu[2:end-1, :]
    @. h.gv *= domainMask.gv[:, 2:end-1]

    # preconditioner and shadow arrays
    e = map(zero, ϕ)
    r̄ = map(zero, r)
    ϕ̄ = map(zero, ϕ)
    ∂ₓϕ̄ = map(zero, ∂ₓϕ)
    ∂ᵧϕ̄ = map(zero, ∂ᵧϕ)
    qx̄ = map(zero, qx)
    qȳ = map(zero, qy)
    D̄ = map(zero, D)
    ρigH̄ = map(zero, ρigH)
    ρwgB̄ = map(zero, ρwgB)
    N̄eff = map(zero, Neff)
    h̄ = map(zero, h)
    m̄ = map(zero, m)
    u_b̄ = map(zero, u_b)
    I = coloring(nx, ny)

    Δt = params.dt*params.sec_per_year;  # time stepsize [s]

    # Convergence tracking
    apt_iters = []    # not a fan of this, it appends with push later on, I'd prefer to define the array size right away (should be possible with max iter)
    apt_errs  = []
    τ         = 0.0   # [s] for monitoring conservation
    Δτ        = 1.0   # this can probably be defined in the struct
    damp      = 0.0

    # monitor conservation
    dA = dx * dy
    ∫hwb_0dA = sum(h.gh) * dA

    # store water sheet thickness from previous time step
    h_prev = deepcopy(h.gh)
    cond = 1.0 .* (h_prev .- h_r .< 1e-10)  # 1 where h_prevc < h_r, else 0

    for iter in 1:maxiter
        # update hydraulic potential solution
        @. ϕ += Δτ * s

        # update effective pressure
        @. Neff = ρigH + ρwgB - ϕ
        @. Neff *= (!floatMask)

        # update water sheet thickness h
        @. h.gh = max((h_prev + Δt * (u_b / l_r) * h_r * cond) / (1 + Δt * (Ã * abs(Neff)^(params.glen_n - 1) * (Neff) + (u_b / l_r) * cond)), 0.0)
        @. h.gh *= domainMask.gh * (!floatMask)
        @. h.gu = abs(h.gh[2:end, :] + h.gh[1:end-1, :]) * 0.5
        @. h.gv = abs(h.gh[:, 2:end] + h.gh[:, 1:end-1]) * 0.5
        @. h.gu *= domainMask.gu[2:end-1, :];
        @. h.gv *= domainMask.gv[:, 2:end-1];

        if do_autotune && (iter % n_autotune == 0)
            @. z0 = z
        end

        # compute residual
        residual!(r, ϕ, ∂ₓϕ, ∂ᵧϕ, qx, qy, D, ρigH, ρwgB, Neff, h, k, α, β, ϕ_fixed, u_b, h_r, l_r, Ã, params.glen_n, dx, dy, m, floatMask, domainMask)

        # apply preconditioner to the residual
        if iter < 10
            # diagonal of J
            for group in I
                fill!(e, 0.0)
                e[group] .= 1.0
                jvp!(e, r̄, r, ϕ, ϕ̄, ∂ₓϕ, ∂ₓϕ̄, ∂ᵧϕ, ∂ᵧϕ̄, qx, qx̄, qy, qȳ, D, D̄, ρigH, ρigH̄, ρwgB, ρwgB̄, Neff, N̄eff, h, h̄, k, α, β, ϕ_fixed, u_b, u_b̄, h_r, l_r, Ã, params.glen_n, dx, dy,
                        m, m̄, floatMask, domainMask)
                @. P[group] = 1.0 / (abs(r̄[group]) + ɛ)
            end
        end

        @. z = P * r

        if do_autotune && (iter % n_autotune == 0)
            # compute A = abs(dot(s, (z - z0))) / dot(s, s)
            A = abs(sum(s .* (z - z0))) / (sum(s .* s) + ɛ)

            damp = 1 + A - 2 * sqrt(A)
        end

        # compute damped search direction
        @. s = damp * s + z

        # monitor convergence
        if iter % ncheck == 0
            err_rel = Δτ * maximum(abs, s) / (maximum(abs, ϕ) + ɛ)
            push!(apt_iters, iter / nx)
            push!(apt_errs, err_rel)

            if err_rel < etol
                # println("GlaDS converged.")
                break
            end
        end

        # advance pseudo time for convergence check
        τ += Δτ
    end

    if update_basal_water_thickness
        gh.basal_water_thickness .= h.gh
    end
    gh.hydraulic_potential_b .= ϕ
    @. Neff *= domainMask.gh * gh.grounded_fraction
    gh.effective_pressure .= max.(basal_hydrology.min_effective_pressure, Neff)

    return model
end