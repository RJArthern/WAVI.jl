"""
    function

"""
function coupled_driver()

### Grid and boundary conditions ###
#nx = 64
#ny = 8
#nx = 160
nx  = 408
ny = 42
nσ = 4
#nσ = 30
x0 = 0.0
y0 = 0.0
dx = 2000.0
dy = 2000.0
dt_s = unset
dt = dt_s/(3600*24*365)
step_thickness=false

h_mask_m=Array{Float64}(undef, nx, ny);
read!("hmask400_40_24_2_edit.bin", h_mask_m)
h_mask=BitArray{2}(undef, nx, ny);
h_mask[findall(x -> x==0,h_mask_m)].=false
h_mask[findall(x -> x==1,h_mask_m)].=true
#Homogenous Dirichlet boundary conditions
#u_iszero=falses(nx+1,ny)
#u_iszero[1,:].=true
#v_iszero=falses(nx,ny+1)
#v_iszero[:,1].=true
#v_iszero[:,end].=true

#orginial
u_iszero=falses(nx+1,ny)
##u_iszero[1,:].=true
u_iszero[9,:].=true
u_iszero[end,:].=true
##u_iszero[end-1,:].=true
#u_iszero[:,1].=true #no for no slip
#u_iszero[:,end].=true #
v_iszero=falses(nx,ny+1)
v_iszero[:,1].=true
v_iszero[:,2].=true
#v_iszero[9,:].=true
v_iszero[:,end].=true
v_iszero[:,end-1].=true
v_iszero[:,end-2].=true
#v_iszero[1,:].=true #added
#v_iszero[end,:].=true #added
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



#Bed 
bed = WAVI.mismip_plus_bed #function definition (hard coded into WAVI for ease, can also read bed here)
#bed_elevation= Array{Float64}(undef, nx, ny);
#read!("bathyREAL_408.box", bed_elevation)
#bed_elevation .= ntoh.(bed_elevation)

#solver parameters
maxiter_picard = 1
solver_params = SolverParams(maxiter_picard = maxiter_picard)

#physical parameters
params = Params() #because accumulation rate commented out, but would set here

#make the model
model = Model(grid = grid,
                    bed_elevation = bed, 
                    params = params, 
                    solver_params = solver_params)


#timestepping parameters
n_iter0 = 0 #CHANGE ME FOR A PICKUP
step_thickness = false #default = true!
dt = 0.1
end_time = 40.
chkpt_freq = 100.
pchkpt_freq = 200.
timestepping_params = TimesteppingParams(n_iter0 = n_iter0, 
                                        dt = dt, 
                                        end_time = end_time, 
                                        chkpt_freq = chkpt_freq, 
                                        pchkpt_freq = pchkpt_freq,
                                        step_thickness = step_thickness)

#### output parameters ###
#outputs = (h = model.gh.h, u = model.gu.u);
#output_freq = 100.
#output_params = OutputParams(outputs = outputs, 
#                        output_freq = output_freq,
#                        format = "mat",
#                        dump_vel = true)
output_params = OutputParams(dump_vel = true) #no outputs, except for dumping the velocity at the end of the run

### set up the simulation ### 
#!! This does a pickup if n_iter0 >0, does not rebuild the model unless flag activated !!
simulation = Simulation(model = model, 
                        timestepping_params = timestepping_params, 
                        output_params = output_params)

### set the thickness if we're doing a pickup
if n_iter0 > 0
    read!(string("Streamice_Thickness_out.data"), thickness_file)
    thickness_file= Array{Float64}(undef, simulation.grid.nx, simulation.grid.ny);
    thickness_file .= ntoh.(thickness_file)
    simulation = @set simulation.gh.h .= thickness_file
end

### run the simulation ###
run_simulation!(simulation)
return nothing

end