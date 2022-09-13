using WAVI 
function driver()

#
#Grid and boundary conditions
#
nx = 160
ny = 20
nσ = 4
x0 = 0.0
y0 = -40000.0
dx = 4000.0
dy = 4000.0
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

#
#Bed 
#
function modified_mismip_plus_bed(x,y)
    xbar = 260000.0
    b0 = -150.0; b2 = -780.0; b4 = 343.91; b6 = -43.6;
    wc = 24000.0; fc = 4000.0; dc = 500.0
    bx(x)=b0+b2*(x/xbar)^2+b4*(x/xbar)^4+b6*(x/xbar)^6
    by(y)= dc*( (1+exp(-2(y-wc)/fc))^(-1) + (1+exp(2(y+wc)/fc))^(-1) )
    b = max(bx(x) + by(y), -720.0)
    return b
end
bed = modified_mismip_plus_bed;


#
#solver parameters
#
maxiter_picard = 1
solver_params = SolverParams(maxiter_picard = maxiter_picard)

#
#Physical parameters
#
accumulation_rate = 0.3
params = Params(accumulation_rate = accumulation_rate)

#
#Initial  Conditionns
# 
h_init = Array{Float64}(undef, nx,ny);
thickness_file = joinpath(dirname(@__FILE__), "ATTR_402_final_thickness.bin")
read!(thickness_file, h_init);
h_init = ntoh.(h_init)
initial_conditions = InitialConditions(initial_thickness = h_init);

#
# Melt rate
#
melt = QuadraticMeltRate(γT = 1.21*1e-3);

#
#make the model
#
model = Model(grid = grid,
              bed_elevation = bed, 
              params = params, 
              solver_params = solver_params,
              initial_conditions = initial_conditions,
              melt_rate = melt);

#
#timestepping parameters
#
niter0 = 0
dt = 0.5
end_time = 100.
chkpt_freq = 1000.
pchkpt_freq = 1000.
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
	   grfrac = model.fields.gh.grounded_fraction,
        m = model.fields.gh.basal_melt)

output_freq = 10.
folder = joinpath(dirname(@__FILE__), "outputs/");
println(folder)
isdir(folder) && rm(folder, force = true, recursive = true)
mkdir(folder) #make a clean folder for outputs
output_params = OutputParams(outputs = outputs, 
                            output_freq = output_freq,
                            output_format = "mat",
                            output_path = folder,
                            dump_vel = false,
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

driver()

