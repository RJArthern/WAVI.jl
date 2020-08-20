"""
Iceberg(n_timesteps)

Driver routine for iceberg simulation using WAVI.
Defines parameters, initialises state, and runs n_timesteps.
"""
function Iceberg_test(n_timesteps=1000)

#Number of timesteps between plots
n_steps_plot=100

# Iceberg
nx = 20
ny = 20
nσ = 4
x0 = -50000.0
y0 = -50000.0
dx = 5000.0
dy = 5000.0
dt = 0.1
xx=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]
yy=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]

#Model domain mask on h-grid
h_mask=sqrt.(xx.^2+yy.^2) .< 45000
starting_thickness=zeros(nx,ny)
starting_thickness[h_mask] .= 200.0
accumulation_rate=0.3



#Homogenous Dirichlet boundary conditions
u_iszero=falses(nx+1,ny)
u_iszero[div(nx,2),:].=true
v_iszero=falses(nx,ny+1)
v_iszero[:,div(ny,2)].=true


bed_elevation=-500.0.*ones(nx,ny)

maxiter_picard=1 #No need for Picard iteration if runninng to steady state

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
            maxiter_picard=maxiter_picard
            )

wavi=start(params)

for i=1:n_timesteps
run!(wavi)
print("\r")
print("Completed: ",round(i*params.dt),
      " years (", round(100.0*i/n_timesteps),"%)  ")
if mod(i,n_steps_plot)==0
  #IJulia.clear_output(true)
  plot_output(wavi)
end
end

return wavi
end
