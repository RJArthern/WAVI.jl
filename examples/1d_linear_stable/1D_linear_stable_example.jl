#Example demonstrating how the sheet/shelf system reaches a stable equilibrium
#on a prograde bed when the ice sheet is one dimensional. This code runs the
#simulation and plots the grounding line as a function of t, demonstrating
#convergence. Note that this is not a simulation of actualy grounding line
#position as only a single picard iteration is included in each time-step
#(running to steady state). Note that this prediction does not agree with Schoof
#2007 boundary layer analysis because we include vertical shearing and our
#effective beta includes variations in z.

function oneD_linear_stable(n_timesteps = 40000)

tic = time()

#parameters
nx = 1000            #number of grid points in x-direction
ny = 3              #number of grid points in y-ndirection
nσ = 4              #number of z-levels
x0 = 0.0            #minimum x-value
x1 = 640.0e3        #maximum x-value
y0 = -40.0e3        #minimum y-value
y1 = 40.e3;         #maximum y-value
dx = (x1 - x0)/nx;  #grid spacing x
dy = (y1 - y0)/ny;  #grid spacing y
dt = 0.5            #timestep
accumulation_rate = 0.3
starting_thickness=200.0.*ones(nx,ny) #initial condition
maxiter_picard = 1 #No need for Picard iteration if runninng to steady state
nsteps_out = 50 #number of steps per output point
glen_reg_strain_rate = 1e-7

#grids
xx=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]
yy=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]
h_mask=trues(nx,ny) #consider ice everywhere in domain

#linear bed elevation
function bed_elev_function(x,y)
bx(x) = 100 + (-382.8 - 100)*x / (3.5*1e5); #linear interpolation of the first 400 or so km of MISMIP+ domain
by(y) = 0
b = bx(x) + by(y)
return b
end
bed_elevation=bed_elev_function.(xx,yy)

#initialize output of time and grounding line
t = zeros((1, floor(Int, n_timesteps/nsteps_out)))
xgl = zeros((1, floor(Int, n_timesteps/nsteps_out)))
ugl = zeros((1, floor(Int, n_timesteps/nsteps_out)))
hgl = zeros((1, floor(Int, n_timesteps/nsteps_out)))

#boundary conditions: Homogenous Dirichlet
u_iszero=falses(nx+1,ny)
u_iszero[1,:].=true
v_iszero=falses(nx,ny+1)
v_iszero[:,1].=true
v_iszero[:,end].=true

#set up params structure
params = Params(nx=nx,
            ny=ny,
            nσ=nσ,
            x0=x0,
            y0=y0,
            dx=dx,
            dy=dy,
            dt=dt,
            bed_elevation=bed_elevation,
            starting_thickness=starting_thickness,
            accumulation_rate=accumulation_rate,
            h_mask=h_mask,
            u_iszero=u_iszero,
            v_iszero=v_iszero,
            maxiter_picard=maxiter_picard,
            glen_reg_strain_rate = glen_reg_strain_rate;
            )
wavi=start(params)

#perform time-stepping
count = 1 #out time-steps counter
for i=1:n_timesteps
      run!(wavi)
      print("\r")
      print("Completed: ",round(i*params.dt),
      " years (", round(100.0*i/n_timesteps),"%)  ")

      #store time and GL position
      if mod(i,nsteps_out)==0
            t[count] = i*params.dt;
            xgl[count] = WAVI.get_GLx(wavi)[2]

            mxval, mxindx= findmax(-broadcast(abs, wavi.gh.haf; ); dims = 1); #returns index of grid point with haf closest to zeros
            u_gl_all = wavi.gh.av_speed[mxindx]
            h_gl_all = wavi.gh.h[mxindx]
            hgl[count] = h_gl_all[2]
            ugl[count] = u_gl_all[2]

            count = count + 1
      end
end

toc = time()
print("Integration took ",round(toc - tic), " secs \n")

return t, xgl, wavi
end

#test run: will not converge as not run for long enough
#t, xgl, wavi, ugl, hgl = oneD_linear_stable(20)

#for convergence, need long time run
t, xgl, wavi = oneD_linear_stable()

#plot outcome
Plots.plot(t', xgl', st = :scatter, legend= false)
