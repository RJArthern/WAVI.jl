using MPI
using WAVI 

function MISMIP_PLUS_GRID(;
        nx = 80,
        ny = 10,
    )
    #Grid and boundary conditions
    nσ = 4
    x0 = 0.0
    y0 = -40000.0
    dx = 8000.0
    dy = 8000.0
    h_mask=trues(nx,ny)
    u_iszero = falses(nx+1,ny); u_iszero[1,:].=true
    v_iszero=falses(nx,ny+1); v_iszero[:,1].=true; v_iszero[:,end].=true
    grid = Grid(nx = nx, 
                ny = ny,   
                nσ = nσ, 
                x0 = x0, 
                y0 = y0, 
                dx = dx, 
                dy = dy,
                h_mask = h_mask, 
                u_iszero = u_iszero, 
                v_iszero = v_iszero)
    return grid
end

function MISMIP_PLUS(;
        folder = "outputs",
        grid = MISMIP_PLUS_GRID(),
        spec = BasicSpec(),
    )
    #Bed 
    bed = WAVI.mismip_plus_bed #function definition

    #solver parameters
    maxiter_picard = 1
    solver_params = SolverParams(maxiter_picard = maxiter_picard)

    #Physical parameters
    default_thickness = 100.0 #set the initial condition this way
    accumulation_rate = 0.3
    default_temperature = 265.700709
    params = Params(default_thickness = default_thickness, 
                    accumulation_rate = accumulation_rate,
                    default_temperature = default_temperature)

    #make the model
    model = Model(grid, bed, spec;
                  params = params, 
                  solver_params = solver_params)

    #timestepping parameters
    niter0 = 0
    dt = .5
    end_time = 400.
    chkpt_freq = 20.
    timestepping_params = TimesteppingParams(niter0 = niter0, 
                                             dt = dt, 
                                             end_time = end_time, 
                                             chkpt_freq = chkpt_freq)

    outputs = (
        h = "global_fields.gh.h",
        u = "global_fields.gh.u",
        v = "global_fields.gh.v",
        b = "global_fields.gh.b",
        grfrac = "global_fields.gh.grounded_fraction",
        mpi_rank = "global_fields.gh.mpi_rank",
    )
    output_freq = 100.
    output_params = OutputParams(outputs,
                                 output_path = folder,
                                 output_freq = output_freq,
                                 output_format = "jld2",
                                 zip_format = "nc",
                                 output_start = true)
    
    simulation = Simulation(model = model, 
                            timestepping_params = timestepping_params,
                            output_params = output_params)
    
    run_simulation!(simulation)
    return simulation
end

if abspath(PROGRAM_FILE) == @__FILE__
    # GPU is opt-in so this script can still run BasicSpec / ThreadedSpec /
    # MPISpec on a GPU node. Example: `WAVI_USE_GPU=1 julia --project=<project_dir> MISMIP_PLUS.jl`.
    # Multiple GPUs: `WAVI_USE_GPU=1 mpiexecjl -n N julia --project=<project_dir> MISMIP_PLUS.jl`.
    use_gpu = get(ENV, "WAVI_USE_GPU", "") == "1"

    # Initialise MPI before CUDA. Some MPI stacks reshuffle visible GPUs if
    # CUDA is already initialised, which breaks one-GPU-per-rank pinning.
    MPI.Init()
    if use_gpu
        try
            using CUDA
        catch
            println("CUDA.jl is not in this project; skipping GPUSpec / MPI+GPU.")
            println("Add CUDA to the active project (not to WAVI [deps]). See docs/src/gpu_setup.md.")
            exit(0)
        end
        if !(gpu_device() isa GPU)
            println("No functional GPU; skipping GPUSpec / MPI+GPU.")
            exit(0)
        end
    end

    if MPI.Comm_size(MPI.COMM_WORLD) > 1
        grid = MISMIP_PLUS_GRID()
        child = use_gpu ? GPU() : CPU()
        mpi_spec = MPISpec(
            MPI.Comm_size(MPI.COMM_WORLD), 1, 2, grid;
            pou = true,
            niterations = 5,
            child_architecture = child,
        )
        folder = use_gpu ? "outputs/mpi_gpu" : "outputs/mpi"
        MISMIP_PLUS(folder = folder, grid = grid, spec = mpi_spec)
    elseif use_gpu
        grid = MISMIP_PLUS_GRID()
        MISMIP_PLUS(folder = "outputs/gpu", grid = grid, spec = GPUSpec())
    elseif Threads.nthreads() > 1
        grid = MISMIP_PLUS_GRID()
        threaded_spec = ThreadedSpec(ngridsx=Threads.nthreads(), ngridsy=1, overlap=2, niterations=5)
        MISMIP_PLUS(
            folder = "outputs/thread",
            grid = grid,
            spec = threaded_spec,
        )
    else
        MISMIP_PLUS(folder = "outputs/serial")
    end
end
