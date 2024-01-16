using WAVI 
function driver()


#
#Grid and boundary conditions
#
nx = 398
ny = 250
nσ = 12
x0 = -1789000.0
y0 = -833000.0
dx = 2000.0
dy = 2000.0

h_mask=Array{Float64}(undef,nx,ny);
read!("Inverse_2km_h_mask_clip_BedmachineV3_FULL_stripe_fix.bin",h_mask)
h_mask.=ntoh.(h_mask)

#
# boundary conditions
#
u_iszero=Array{Float64}(undef,nx+1,ny);
read!("Inverse_2km_uiszero_clip_BedmachineV3_FULL_stripe_fix.bin",u_iszero)
u_iszero.=ntoh.(u_iszero)

v_iszero=Array{Float64}(undef,nx,ny+1);
read!("Inverse_2km_viszero_clip_BedmachineV3_FULL_stripe_fix.bin",v_iszero)
v_iszero.=ntoh.(v_iszero)

sigma_grid=Array{Float64}(undef,nσ);
read!("Inverse_2km_sigma_grid_BedmachineV3_FULL_stripe_fix.bin",sigma_grid)
sigma_grid.=ntoh.(sigma_grid)

grid = Grid(nx = nx,
            ny = ny,
            nσ = nσ,
            x0 = x0,
            y0 = y0,
            dx = dx,
            dy = dy,
            h_mask = h_mask,
            u_iszero = u_iszero,
            v_iszero = v_iszero,
            σ = sigma_grid)

#
#input files
#
bed=Array{Float64}(undef,nx,ny);
read!("Inverse_2km_bed_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",bed)
bed.=ntoh.(bed)

h=Array{Float64}(undef,nx,ny);
read!("Inverse_2km_thickness_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",h)
h.=ntoh.(h)

viscosity=Array{Float64}(undef,nx,ny,nσ);
read!("Inverse_2km_viscosity3D_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",viscosity)
viscosity.=ntoh.(viscosity)

temp=Array{Float64}(undef,nx,ny,nσ);
read!("Inverse_2km_3Dtemp_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",temp)
temp.=ntoh.(temp)

damage=Array{Float64}(undef,nx,ny,nσ);
read!("Inverse_2km_damage3D_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",damage)
damage.=ntoh.(damage)

weertman_c=Array{Float64}(undef,nx,ny);
read!("Inverse_2km_WeertmanC_clip_adjusted_noNan_BedmachineV3_FULL_stripe_fix.bin",weertman_c)
weertman_c.=ntoh.(weertman_c)

accumulation_rate=Array{Float64}(undef,nx,ny);
read!("Inverse_2km_accumulation_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",accumulation_rate)
accumulation_rate.=ntoh.(accumulation_rate)

dhdt=Array{Float64}(undef,nx,ny);
read!("Inverse_2km_dhdt_clip_noNan_BedmachineV3_FULL_stripe_fix.bin",dhdt)
dhdt.=ntoh.(dhdt)

initial_conditions = InitialConditions(initial_thickness = h,
                                        initial_viscosity = viscosity,
                                        initial_temperature = temp,
                                        initial_damage = damage)



#
#solver parameters
#
maxiter_picard = 1
#parallel_spec = SharedMemorySpec(ngridsx=2,ngridsy=1,overlap=1,damping=0.0,niterations=1)
parallel_spec = BasicParallelSpec()
tol_picard = 1.0e-4;
solver_params = SolverParams(maxiter_picard = maxiter_picard, tol_picard = tol_picard)

params = Params(accumulation_rate = accumulation_rate,
				  weertman_c = weertman_c)

bump_amplitude = 0.
bump_width = 2.5
bump_time = 245. #put beyond end of simulation
trend_onset = 260. #beyond end of simulation
pc_max = -399.0
pc_min = -401.0 #warm conditions
pc_max = -599.0
pc_min = -601.0 #warm conditions
pw     = 400.0
rf_threshold = 0.01 #collpase everything

idealized_anthro_melt_rate = IdealizedAnthroMeltRate(bump_amplitude = bump_amplitude,
bump_width = bump_width,
bump_time = bump_time,
per_century_trend = 0.0,
trend_onset = trend_onset,
pc_max = pc_max,
pc_min = pc_min,
M = 2.6,
random_seed = 1,
rf_threshold = rf_threshold,
pw = pw)


model = Model(grid = grid,
            bed_elevation = bed,
            params = params,
            solver_params = solver_params,
            initial_conditions= initial_conditions, 
            parallel_spec = parallel_spec,
	    melt_rate = idealized_anthro_melt_rate)


#
#timestepping parameters
#
niter0 = 0
dt = 0.05
end_time = dt
#end_time = 1000.
chkpt_freq = 10.0
pchkpt_freq = 10.0
timestepping_params = TimesteppingParams(niter0 = niter0, 
                                         dt = dt, 
                                         end_time = end_time, 
                                         chkpt_freq = chkpt_freq, 
                                          pchkpt_freq = pchkpt_freq)

#
#output parameters
#
outputs = (h = model.fields.gh.h,
           u = model.fields.gh.u,
           v = model.fields.gh.v,
           b = model.fields.gh.b,
           s = model.fields.gh.s,
	       a = model.fields.gh.accumulation,
           grfrac = model.fields.gh.grounded_fraction,
           m = model.fields.gh.basal_melt)

output_freq = 5.0
output_params = OutputParams(outputs = outputs, 
                            output_freq = output_freq,
                            output_format = "mat",
                            zip_format = "nc",
                            output_start = true)

#
# assemble the simulation
#
simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params, 
                        output_params = output_params)
                
#
#perform the simulation
#
run_simulation!(simulation)

return simulation
end

simulation = driver()

