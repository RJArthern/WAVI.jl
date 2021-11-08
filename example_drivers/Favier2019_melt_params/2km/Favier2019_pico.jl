using WAVI, Plots
# Grid and boundary conditions
nx = 320
ny = 40
nσ = 4
x0 = 0.0
y0 = -40000.0
dx = 2000.0
dy = 2000.0
h_mask = trues(nx, ny)
u_iszero = falses(nx + 1, ny); u_iszero[1,:] .= true
v_iszero = falses(nx, ny + 1); v_iszero[:,1] .= true; v_iszero[:,end] .= true
grid = Grid(nx=nx, 
            ny=ny,   
            nσ=nσ, 
            x0=x0, 
            y0=y0, 
            dx=dx, 
            dy=dy,
            h_mask=h_mask, 
            u_iszero=u_iszero, 
            v_iszero=v_iszero)

# Bed 
bed = WAVI.mismip_plus_bed # function definition

# Inputing thickness profile, reading from binary file 
#fname = "examples\\Favier2019_melt_params\\data\\MISMIP_ice0_2km_SteadyThickness.bin";
fname = joinpath(dirname(@__FILE__), "data",  "MISMIP_ice0_2km_SteadyThickness.bin")
h = Array{Float64,2}(undef, nx, ny)
read!(fname, h)
h = ntoh.(h)

# make the model
if ~(@isdefined model)
initial_conditions = InitialConditions(initial_thickness=h) # set thickness 
model = Model(grid=grid,
            bed_elevation=bed, 
            initial_conditions=initial_conditions,
            solver_params=SolverParams(maxiter_picard=1) # we'll update state to get basal melt, but don't need other info accurate
            ) 

# update the configuration
end
update_state!(model)

ice_front_mask = zeros(nx,ny);
ice_front_mask[end,:] .= 1;
pico = PICO(ice_front_mask = ice_front_mask, 
            T0 =1.2, 
            S0 = 34.6, 
            use_box_mean_depth = true, 
            γT = 3.e-5);

#return the gl mask
ice_front_mask = zeros(nx,ny);
ice_front_mask[end,:] .= 1;
zb = model.fields.gh.b .* (model.fields.gh.grounded_fraction .== 1) + - pico.ρi / pico.ρw .* model.fields.gh.h .* (model.fields.gh.grounded_fraction .< 1)
gl_mask = WAVI.get_gl_mask(model.fields.gh.grounded_fraction)
boxes = WAVI.get_pico_boxes(gl_mask,ice_front_mask, model.fields.gh.grounded_fraction, pico.nbox )

WAVI.update_melt_rate!(pico, model.fields,model.grid)

#contour plot the boxes
x = model.grid.xxh[:,1]; y = model.grid.yyh[1,:];
plt1 = contour(x ./ 1e3,y / 1e3, boxes', 
                fill=true, 
                linewidth=0, 
                colorbar=true,
                colorbar_title="pico boxes",
                framestyle=:box)
xlims!((420, 640))
xlabel!("x (km)")
ylabel!("y (km)")
plt1 

#contour plot melt rate
m = deepcopy(model.fields.gh.basal_melt);
m[model.fields.gh.grounded_fraction .== 1.] .= NaN;
msat = deepcopy(m)
msat[msat .> 50] .= 50
#msat[zb .> -300] .= NaN;
plt2 = contour(x ./ 1e3,y / 1e3, msat', 
                fill=true, 
                linewidth=0, 
                colorbar=true,
                colorbar_title="melt rate (m/yr)",
                framestyle=:box, 
                colorscale = "Jet")
xlims!((420, 640))
xlabel!("x (km)")
ylabel!("y (km)")

m = deepcopy(model.fields.gh.basal_melt);
idx = (zb .< -300) .& (model.fields.gh.basal_melt .> 0)
mean_melt = sum(m[idx])./length(m[idx])
println("mean melt rate for shelf points below -300m depth is $mean_melt")