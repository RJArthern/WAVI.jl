#This will be my relax

using WAVI 
using Printf
using ImageFiltering

function driver()

#
#Grid and boundary conditions
#
nx = 164
ny = 192
nσ = 12
x0 = -1802500.0
y0 = -847500.
dx = 5000.0
dy = 5000.0

#h_mask=trues(nx,ny)
#u_iszero = falses(nx+1,ny); u_iszero[1,:].=true
#v_iszero=falses(nx,ny+1); v_iszero[:,1].=true; v_iszero[:,end].=true

 #println("loading h_mask")
h_mask=Array{Float64}(undef,nx,ny);
read!("..//..//..//..//..//data//PIGTHW5km//Inverse_5km_h_mask_clip.bin",h_mask)
#read!("//data//PIGTHW5km//Inverse_5km_h_mask_clip.bin",h_mask)
h_mask.=ntoh.(h_mask)

 #println("loaded h_mask")

u_iszero=Array{Float64}(undef,nx+1,ny);
read!("..//..//..//..//..//data//PIGTHW5km//Inverse_5km_uiszero_clip.bin",u_iszero)
u_iszero.=ntoh.(u_iszero)

v_iszero=Array{Float64}(undef,nx,ny+1);
read!("..//..//..//..//..//data//PIGTHW5km//Inverse_5km_viszero_clip.bin",v_iszero)
v_iszero.=ntoh.(v_iszero)

sigma_grid=Array{Float64}(undef,nσ);
read!("..//..//..//..//..//data//PIGTHW5km//Inverse_5km_sigma_grid.bin",sigma_grid)
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
#Bed 
#
bed=Array{Float64}(undef,nx,ny);
read!("..//..//..//..//..//data//PIGTHW5km//Inverse_5km_bed_clip_noNan.bin",bed)
bed.=ntoh.(bed)
#bed = WAVI.mismip_plus_bed #function definition

h=Array{Float64}(undef,nx,ny);
read!("..//..//..//..//..//data//PIGTHW5km//Inverse_5km_thickness_clip_noNan.bin",h)
h.=ntoh.(h)

viscosity=Array{Float64}(undef,nx,ny,nσ);
read!("..//..//..//..//..//data//PIGTHW5km//Inverse_5km_viscosity3D_clip_noNan.bin",viscosity)
viscosity.=ntoh.(viscosity)

#
temp=Array{Float64}(undef,nx,ny,nσ);
read!("..//..//..//..//..//data//PIGTHW5km//Inverse_5km_3Dtemp_clip.bin",temp)
temp.=ntoh.(temp)

damage=Array{Float64}(undef,nx,ny,nσ);
read!("..//..//..//..//..//data//PIGTHW5km//Inverse_5km_damage3D_clip_noNan.bin",damage)
damage.=ntoh.(damage)

 
weertman_c=Array{Float64}(undef,nx,ny);
read!("..//..//..//..//..//data//PIGTHW5km//Inverse_5km_WeertmanC_clip_adjusted_noNan.bin",weertman_c)
weertman_c.=ntoh.(weertman_c)

accumulation_rate=Array{Float64}(undef,nx,ny);
read!("..//..//..//..//..//data//PIGTHW5km//Inverse_5km_accumulation_clip_noNan.bin",accumulation_rate)
accumulation_rate.=ntoh.(accumulation_rate)

 
dhdt=Array{Float64}(undef,nx,ny);
read!("..//..//..//..//..//data//PIGTHW5km//Inverse_5km_dhdt_clip_noNan.bin",dhdt)
dhdt.=ntoh.(dhdt)

initial_conditions = InitialConditions(initial_thickness = h,
                                        initial_viscosity = viscosity,
                                        initial_temperature = temp,
                                        initial_damage = damage)



#solver parameters
#
maxiter_picard = 1
tol_picard = 1.0e-4
solver_params = SolverParams(maxiter_picard = maxiter_picard,
                            tol_picard = tol_picard)

#
#Physical parameters
#default_thickness = 100.0 #set the initial condition this way
#accumulation_rate = 0.3


#Because this is a relaxation, we do not evolve the shelves:
evolveShelves = false

#w_centered=centered([1 1 1; 1 1 1; 1 1 1])

#In relaxation, we modfiy the accumulation and run it with a_rlx=a-dh/dt:

accumulation_rate = accumulation_rate-dhdt 

params = Params(accumulation_rate = accumulation_rate,
                weertman_c = weertman_c,
                evolveShelves = evolveShelves)


 #make the model
#
@printf "Starting to make the model"
model = Model(grid = grid,
              bed_elevation = bed, 
              params = params, 
              solver_params = solver_params,
              initial_conditions= initial_conditions) 

#
@sprintf "The model is made"
#timestepping parameters
niter0 = 0
dt = 0.025
end_time = 1.0
chkpt_freq = 1.0
pchkpt_freq = 10.0
timestepping_params = TimesteppingParams(niter0 = niter0, 
                                           dt = dt, 
                                           end_time = end_time, 
                                           chkpt_freq = chkpt_freq, 
                                           pchkpt_freq = pchkpt_freq)

##output parameters
folder = "outputs_5km_relax_test"
isdir(folder) && rm(folder, force = true, recursive = true)
mkdir(folder) #make a clean folder for outputs
outputs = (h = model.fields.gh.h,
            u = model.fields.gh.u,
            v = model.fields.gh.v,
            b = model.fields.gh.b,
            h_mask = model.fields.gh.mask,
            visc_av = model.fields.gh.ηav,
            weertman_c = model.fields.gh.weertman_c,
            grounded_frac = model.fields.gh.grounded_fraction,
            av_speed = model.fields.gh.av_speed,
            melt = model.fields.gh.basal_melt,
            dhdt = model.fields.gh.dhdt,
            beta2 = model.fields.gh.β,
            accumulation = model.fields.gh.accumulation,
            bed_speed = model.fields.gh.bed_speed,
            tau_bed = model.fields.gh.τbed)
             #output velocities and thickness
output_freq = 0.1
output_params = OutputParams(outputs = outputs, 
                           output_path = folder,
                           output_freq = output_freq,
                           output_format = "mat",
                           zip_format = "nc")
   
                           
@printf "About to make simulation"
simulation = Simulation(model = model, 
                      timestepping_params = timestepping_params,
                      output_params = output_params)
           
 @printf "The simulation is made"
#perform the simulation

run_simulation!(simulation)

@printf "The simulation has been run"

return simulation


end

