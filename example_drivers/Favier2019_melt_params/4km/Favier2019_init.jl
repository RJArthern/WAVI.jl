using WAVI, Plots
"""
Produce a plot of the melt rate in MISMIP for specified melt rate parametrization (c.f. figure 4 in Favier 2019 10.5194/gmd-12-2255-2019)
***Options***
    - PICO_nboxP_zQ    : PICO with P boxes, ambient temperature taken at Qm depth for P = 10, 8, 5, 2 and Q = 500, 700
                        (e.g. "PICO_nbox2_z700")
    - PME              : Plume model emulator (using Lazeroms 2018 algorithm for grounding line and basal slope -- doi:10.5194/tc-12-49-2018)
    - MISMIP_1r        : Melt rate according to the MISMIP+ ice 1r experiment (doi: 10.5194/tc-14-2283-2020) scaled to match the mean melt rate
    - QuadL            : Quadratic formulation of melting with local dependency on thermal driving
    - QuadNL           : Quadratic formulation of melting with non-local dependency on thermal driving
"""

melt_rate_model = "QuadL"

function Favier2019_4km_init(melt_model)
# Grid and boundary conditions
nx = 160
ny = 20
nσ = 4
x0 = 0.0
y0 = -40000.0
dx = 4000.0
dy = 4000.0
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
#fname = joinpath(dirname(@__FILE__), "data",  "WAVI_ice0_4km_thick_interpolated.bin")
fname = joinpath(dirname(@__FILE__), "data",  "WAVI_ice0_4km_thickness.bin")
h = Array{Float64,2}(undef, nx, ny)
read!(fname, h)
h = ntoh.(h)

# make the model
initial_conditions = InitialConditions(initial_thickness=h) # set thickness 
model = Model(grid=grid,
            bed_elevation=bed,
            melt_rate = melt_model, 
            initial_conditions=initial_conditions,
            solver_params=SolverParams(maxiter_picard=1) 
            ) ;

# update the configuration
update_state!(model)

#contour plot melt rate
m = deepcopy(model.fields.gh.basal_melt)
m[model.fields.gh.grounded_fraction .== 1.] .= NaN
msat = deepcopy(m)
msat[msat .> 50] .= 50
#msat[zb .> -300] .= NaN;
x =  model.grid.xxh[:,1]; y = model.grid.yyh[1,:];
plt = heatmap(x ./ 1e3,y / 1e3, msat', 
                fill=true, 
                linewidth=0, 
                colorbar=true,
                colorbar_title="melt rate (m/yr)",
                framestyle=:box)

xlims!((420, 640))
xlabel!("x (km)")
ylabel!("y (km)")

return model, plt
end

#Need these for PICO 
nx = 160; ny = 20;
ice_front_mask = zeros(nx,ny);
ice_front_mask[end,:] .= 1;

if melt_rate_model == "PICO_nbox10_z700"
    melt_model = PICO(ice_front_mask = ice_front_mask, 
                T0 =1.2, 
                S0 = 34.6, 
                use_box_mean_depth = true, 
                γT = 0.94e-5, 
                nbox = 10);

elseif melt_rate_model == "PICO_nbox8_z700"
    melt_model = PICO(ice_front_mask = ice_front_mask, 
                T0 =1.2, 
                S0 = 34.6, 
                use_box_mean_depth = true, 
                γT = 0.91e-5, 
                nbox = 8);

elseif melt_rate_model == "PICO_nbox5_z700"
    melt_model = PICO(ice_front_mask = ice_front_mask, 
                T0 = 1.2, 
                S0 = 34.6, 
                use_box_mean_depth = true, 
                γT = 0.87e-5, 
                nbox = 5);

elseif melt_rate_model == "PICO_nbox2_z700"
    melt_model = PICO(ice_front_mask = ice_front_mask, 
                T0 = 1.2, 
                S0 = 34.6, 
                use_box_mean_depth = true, 
                γT = 0.85e-5, 
                nbox = 2);               

elseif melt_rate_model == "PME"  
    melt_model = PlumeEmulator(α=1.49)

elseif melt_rate_model == "MISMIP_1r"
    melt_model = MISMIPMeltRateOne(α = 0.184)

elseif melt_rate_model == "QuadL"
    melt_model = QuadraticMeltRate(γT = 0.745*1e-3)
    #melt_model = QuadraticMeltRate(γT = 0.745*1e-3, melt_partial_cell = true)
    
elseif melt_rate_model == "QuadNL"
    melt_model = QuadraticMeltRate(γT = 0.85*1e-3, flocal = false)
    #melt_model = QuadraticMeltRate(γT = 0.85*1e-3, flocal = false, melt_partial_cell = true)

else
    throw(ArgumentError("Specified melt rate model not found"))
end



model, plt = Favier2019_4km_init(melt_model);

m = deepcopy(model.fields.gh.basal_melt);
zb = model.fields.gh.b .* (model.fields.gh.grounded_fraction .== 1) + -  918.0 / 1028.0  .* model.fields.gh.h .* (model.fields.gh.grounded_fraction .< 1)
idx = (model.fields.gh.grounded_fraction .== 0)
mean_melt = sum(m[idx])./length(m[idx])
println("mean melt rate for shelf points is $mean_melt")