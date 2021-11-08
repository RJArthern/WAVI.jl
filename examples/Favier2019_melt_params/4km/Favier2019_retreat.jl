using WAVI, Plots
"""
Run retreat experiment (100 years from MISMIP steady state) 
***Options***
    - PICO_nboxP_zQ    : PICO with P boxes, ambient temperature taken at Qm depth for P = 10, 8, 5, 2 and Q = 500, 700
                        (e.g. "PICO_nbox2_z700")
    - PME              : Plume model emulator (using Lazeroms 2018 algorithm for grounding line and basal slope -- doi:10.5194/tc-12-49-2018)
    - MISMIP_1r        : Melt rate according to the MISMIP+ ice 1r experiment (doi: 10.5194/tc-14-2283-2020) scaled to match the mean melt rate
"""

melt_rate_model = "PME"

function Favier2019_4km_retreat(melt_model)
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
            solver_params=SolverParams(maxiter_picard=10) ,
            params=Params(accumulation_rate=0.3)
            ) ;

#timestepping parameters
niter0 = 0
dt = 0.02
end_time = 10.
chkpt_freq = 0.5
pchkpt_freq = 0.5
timestepping_params = TimesteppingParams(niter0 = niter0, 
                                        dt = dt, 
                                        end_time = end_time, 
                                        chkpt_freq = chkpt_freq, 
                                        pchkpt_freq = pchkpt_freq)

#output parameters
folder =  joinpath(dirname(@__FILE__),"retreat_output")
isdir(folder) && rm(folder, force = true, recursive = true)
mkdir(folder) #make a clean folder for outputs
outputs = (h   = model.fields.gh.h,
            u  = model.fields.gh.u,
            v  = model.fields.gh.v,
            melt = model.fields.gh.basal_melt,
            grfrac = model.fields.gh.grounded_fraction)
output_freq = 0.5
output_params = OutputParams(outputs = outputs, 
                        output_freq = output_freq,
                        output_format = "jld2",
                        output_path = folder,
                        zip_format = "nc")

simulation = Simulation(model = model, 
                    timestepping_params = timestepping_params, 
                    output_params = output_params);
        
#perform the simulation
run_simulation!(simulation);

return simulation
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
    melt_model = MISMIPMeltRateOne(α = 0.32)
else
    throw(ArgumentError("Specified melt rate model not found"))
end



@time simulation = Favier2019_4km_retreat(melt_model);

#m = deepcopy(model.fields.gh.basal_melt);
#zb = model.fields.gh.b .* (model.fields.gh.grounded_fraction .== 1) + -  918.0 / 1028.0  .* model.fields.gh.h .* (model.fields.gh.grounded_fraction .< 1)
#idx = (model.fields.gh.grounded_fraction .== 0)
#mean_melt = sum(m[idx])./length(m[idx])
#println("mean melt rate for shelf points is $mean_melt")