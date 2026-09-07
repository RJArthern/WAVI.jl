using Test
using WAVI
using SparseArrays
using WAVI.Architectures: architecture, zeros_on, array_type

@testset "GPUSpec" begin
    @testset "CPU helpers" begin
        arch = CPU()
        z = zeros_on(arch, Float64, 2, 3)
        @test z isa Array{Float64,2}
        @test size(z) == (2, 3)
        @test all(iszero, z)
        @test architecture(arch) === arch
        @test architecture(BasicSpec()) isa CPU
        @test assign_local_device!(CPU(), 0) === nothing

        grid = Grid(nx = 8, ny = 8, dx = 1.0e3, dy = 1.0e3)
        model = Model(grid = grid, bed_elevation = zeros(grid.nx, grid.ny), verbose = false)
        update_state!(model)
        @test model.fields.stencil_scratch[] !== nothing
        restored = WAVI.Models.restore_pickup_architecture!(model)
        @test restored.fields.stencil_scratch[] === nothing
        update_state!(restored)
        @test all(isfinite, restored.fields.gu.u)

        @testset "update_state stencils match Kronecker" begin
            grid = Grid(nx = 8, ny = 8, dx = 1.0e3, dy = 1.0e3)
            model = Model(grid = grid, bed_elevation = zeros(grid.nx, grid.ny), verbose = false)
            WAVI.Processes.update_surface_elevation!(model)
            gh, gu, gv = model.fields.gh, model.fields.gu, model.fields.gv
            onesvec = ones(length(gh.h))
            denu = gu.samp * (gu.centᵀ * (gh.crop * onesvec))
            denv = gv.samp * (gv.centᵀ * (gh.crop * onesvec))
            h_u = (gu.samp * (gu.centᵀ * (gh.crop * vec(gh.h)))) ./ denu
            h_v = (gv.samp * (gv.centᵀ * (gh.crop * vec(gh.h)))) ./ denv
            WAVI.Processes.update_geometry_on_uv_grids!(model)
            @test gu.h[gu.mask] ≈ h_u
            @test gv.h[gv.mask] ≈ h_v

            gu.u[gu.mask] .= 50.0
            gv.v[gv.mask] .= -30.0

            u_h = reshape(gu.cent * vec(gu.u), size(gh.u))
            v_h = reshape(gv.cent * vec(gv.v), size(gh.v))
            WAVI.Processes.update_velocities_on_h_grid!(model)
            @test gh.u ≈ u_h
            @test gh.v ≈ v_h

            us_u = reshape(gu.crop * (gu.centᵀ * (gh.crop * vec(gh.us))), size(gu.us))
            vs_v = reshape(gv.crop * (gv.centᵀ * (gh.crop * vec(gh.vs))), size(gv.vs))
            WAVI.Processes.update_surface_velocities_on_uv_grid!(model)
            @test gu.us ≈ us_u
            @test gv.vs ≈ vs_v

            dhdt_packed = gh.samp * (
                vec(gh.accumulation) .- vec(gh.basal_melt) .- (
                    (gu.∂x * (gu.crop * (vec(gu.h) .* vec(gu.u)))) .+
                    (gv.∂y * (gv.crop * (vec(gv.h) .* vec(gv.v))))
                )
            )
            WAVI.Processes.update_dhdt!(model)
            @test gh.dhdt[gh.mask] ≈ dhdt_packed
        end
    end

    cuda_ext = Base.get_extension(WAVI, :WAVICUDAExt)
    if cuda_ext === nothing
        @testset "GPUSpec needs a GPU backend" begin
            @test_throws Exception GPUSpec()
        end
    else
        @testset "GPUSpec with CUDA extension" begin
            detected = gpu_device()
            if !(detected isa GPU)
                @test_throws Exception GPUSpec()
            else
                spec = GPUSpec()
                @test spec isa GPUSpec
                @test architecture(spec) isa GPU
                @test array_type(architecture(spec)) !== Array

                grid = Grid(nx = 8, ny = 8, dx = 1.0e3, dy = 1.0e3)
                bed = zeros(grid.nx, grid.ny)
                model = Model(grid = grid, bed_elevation = bed, spec = spec)
                @test architecture(model) isa GPU
                AT = array_type(architecture(spec))
                @test model.fields.gh.h isa AT
                @test model.fields.gu.u isa AT
                @test model.fields.wu.wavelets isa AT
                @test model.fields.gh.samp isa SparseMatrixCSC

                s = WAVI.Utilities.stencil_scratch!(model)
                @test s.haar_u isa AT
                @test s.rhs isa AT

                update_state!(model)
                u_gpu = Array(model.fields.gu.u)
                v_gpu = Array(model.fields.gv.v)

                cpu_model = Model(grid = grid, bed_elevation = bed)
                update_state!(cpu_model)
                @test u_gpu ≈ cpu_model.fields.gu.u rtol = 1e-8 atol = 1e-8
                @test v_gpu ≈ cpu_model.fields.gv.v rtol = 1e-8 atol = 1e-8
                @test Array(model.fields.gh.dhdt) ≈ cpu_model.fields.gh.dhdt rtol = 1e-8 atol = 1e-8
                @test Array(model.fields.gu.us) ≈ cpu_model.fields.gu.us rtol = 1e-8 atol = 1e-8

                @testset "GPUSpec checkpoint pickup and resume" begin
                    dt = 0.1
                    t_mid = 0.2
                    t_end = 0.4
                    mktempdir() do dir
                        function gpu_chkpt_model()
                            return Model(
                                grid = grid,
                                bed_elevation = bed,
                                spec = GPUSpec(),
                                params = Params(accumulation_rate = 0.1),
                                solver_params = SolverParams(maxiter_picard = 1),
                                verbose = false,
                            )
                        end
                        sim_write = Simulation(
                            model = gpu_chkpt_model(),
                            timestepping_params = TimesteppingParams(
                                dt = dt,
                                end_time = t_mid,
                                chkpt_freq = t_mid,
                                chkpt_path = dir,
                            ),
                            output_params = OutputParams(output_path = dir),
                        )
                        run_simulation!(sim_write)
                        n_iter = sim_write.clock.n_iter
                        h0 = Array(sim_write.model.fields.gh.h)
                        u0 = Array(sim_write.model.fields.gu.u)
                        v0 = Array(sim_write.model.fields.gv.v)

                        sim_pick = Simulation(
                            model = gpu_chkpt_model(),
                            timestepping_params = TimesteppingParams(
                                dt = dt,
                                end_time = t_end,
                                niter0 = n_iter,
                                chkpt_path = dir,
                            ),
                            output_params = OutputParams(output_path = dir),
                        )
                        @test sim_pick.clock.n_iter == n_iter
                        @test sim_pick.model.fields.gh.h isa AT
                        @test sim_pick.model.fields.stencil_scratch[] === nothing
                        @test Array(sim_pick.model.fields.gh.h) ≈ h0
                        @test Array(sim_pick.model.fields.gu.u) ≈ u0
                        @test Array(sim_pick.model.fields.gv.v) ≈ v0

                        run_simulation!(sim_pick)
                        @test sim_pick.clock.time ≈ t_end
                        @test sim_pick.model.fields.gh.h isa AT

                        sim_ctrl = Simulation(
                            model = gpu_chkpt_model(),
                            timestepping_params = TimesteppingParams(
                                dt = dt,
                                end_time = t_end,
                                chkpt_path = dir,
                            ),
                            output_params = OutputParams(output_path = dir),
                        )
                        run_simulation!(sim_ctrl)
                        @test Array(sim_pick.model.fields.gh.h) ≈ Array(sim_ctrl.model.fields.gh.h)
                        @test Array(sim_pick.model.fields.gu.u) ≈ Array(sim_ctrl.model.fields.gu.u)
                        @test Array(sim_pick.model.fields.gv.v) ≈ Array(sim_ctrl.model.fields.gv.v)
                    end
                end
            end
        end
    end
end
