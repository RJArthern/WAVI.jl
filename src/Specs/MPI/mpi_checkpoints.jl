

# MPI Checkpointing. Overloads checkpointing routines so that each rank outputs it's own checkpoint file.

 checkpoint_filename(spec::MPISpec,n_iter::Integer) = string("Chkpt_", lpad(n_iter, 10, "0"),"_Rank",lpad(spec.rank, 10, "0"), ".jld2")


 function load_checkpoint(spec::MPISpec,timestepping_params::TimesteppingParams, output_params::OutputParams)

    @unpack comm = spec

    path = checkpoint_path(timestepping_params, output_params)

    n_iter = timestepping_params.niter0
    fname = joinpath(path, checkpoint_filename(spec,n_iter))


    MPI.Barrier(comm)
    dict = load(fname)
    MPI.Barrier(comm)

    if haskey(dict, "simulation")
        # legacy format (pre #105): @save fname simulation saved everything
        # in one bundle, need to unpack it.
        sim_load = dict["simulation"]
        model = sim_load.model
        clock = sim_load.clock
    else
        # new format (post #105): @save fname model=model clock=clock saved separately
        model = dict["model"]
        clock = dict["clock"]
    end
    if clock.n_iter != n_iter
        @warn "Checkpoint iteration $(clock.n_iter) does not match niter0=$n_iter ($(fname))"
    end

    # Check consistency of current spec relative to pickup
    (spec.px == model.spec.px &&
    spec.py == model.spec.py &&
    spec.halo == model.spec.halo &&
    spec.global_size == model.spec.global_size &&
    spec.rank == model.spec.rank &&
    spec.coords == model.spec.coords &&
    spec.top == model.spec.top &&
    spec.right == model.spec.right &&
    spec.bottom == model.spec.bottom &&
    spec.left == model.spec.left) || error("Model spec from pickup is inconsistent with current MPISpec")

    saved_global_fields = model.spec.global_fields

    #Rebuild model using current MPISpec rather than saved spec.
    model = Model(model.grid,
                  model.fields,
                  model.params,
                  model.solver_params,
                  spec,
                  model.shelf_melt_rate,
                  model.surface_mass_balance,
                  model.fracture,
                  model.sliding_law,
                  model.basal_hydrology,
                  model.thermo_dynamics,
                  model.verbose)
    mpi_restore_global_fields!(spec, saved_global_fields, model)
    return model, clock
end

function write_checkpoint!(spec::MPISpec,model, timestepping_params::TimesteppingParams, output_params::OutputParams, clock::Clock)
    path = checkpoint_path(timestepping_params, output_params)

    @unpack rank,comm = spec

    MPI.Barrier(comm)
    if rank==0 && !isdir(path)
        mkpath(path)
    end
    MPI.Barrier(comm)
    
    fname = joinpath(path, checkpoint_filename(spec,clock.n_iter))

    MPI.Barrier(comm)
    with_cleared_stencil_scratch(model) do
        @save fname model=model timestepping_params=timestepping_params clock=clock
    end
    MPI.Barrier(comm)

    @info "Permanent checkpoint at timestep number $(clock.n_iter) — $(fname)"
    return nothing
end