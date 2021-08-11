nx = 320
ny = 40
nσ = 4
x0 = 0.0
y0 = -40000.0
dx = 2000.0
dy = 2000.0
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

#Bed 
bed = WAVI.mismip_plus_bed #function definition

#Inputing thickness profile, reading from binary file 
fname = "examples\\Favier2019_melt_params\\data\\MISMIP_ice0_2km_SteadyThickness.bin";
#fname = joinpath(dirname(@__FILE__), "data",  "MISMIP_ice0_2km_SteadyThickness.bin")
h = Array{Float64, 2}(undef, nx, ny)
read!(fname,h)
h = ntoh.(h)        

#make the model
initial_conditions = InitialConditions(initial_thickness = h) #set thickness 
model = Model(grid = grid,
            bed_elevation = bed, 
            initial_conditions = initial_conditions) #don't care about params or solver params

#embed the ice model with melt rate model
pme = PlumeEmulator(α = 1.46) #change prefactor

#update grounded_fraction
WAVI.update_height_above_floatation!(model)
WAVI.update_grounded_fraction_on_huv_grids!(model) 


melt = zeros(nx,ny);
h = model.fields.gh.h;
grounded_fraction = model.fields.gh.grounded_fraction;
bathy = model.fields.gh.b;
ρi = model.params.density_ice;  
ρw = model.params.density_ocean;
dx = model.grid.dx;
dy = model.grid.dy;

#first steps in set_plume_emulator_melt_rate:
nx, ny = size(h);
zbf = @. -(ρi/ρw)*h.*(1-grounded_fraction) .+ grounded_fraction*bathy; #ice draft if floating everywhere
∂zb∂x, ∂zb∂y = WAVI.get_slope(zbf,dx, dy); #returns the partial derivates of base in both directions

i = 230; j = 20; #pt in shelf
k = 13; 
s = zeros(1, nd); #basal slopes


dirs = WAVI.get_plume_directions(); #array of directions
nd = size(dirs)[1]; #number of search directions
norm_dirs = dirs ./ norm.(eachrow(dirs)); #normalized directions    

#slope 
s[k] = -dot(norm_dirs[k,:], [∂zb∂x[i,j] ∂zb∂y[i,j]]); #note negative sign

