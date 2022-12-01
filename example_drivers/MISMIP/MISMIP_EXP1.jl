using WAVI, Plots 
"""
    MISMIP Experiments EXP1_1: advance to a steady state on a linear slope

"""
function MISMIP_EXP1_1(A)
    #grid and bc
    nx = 300
    ny = 2
    nσ = 4
    x0 = 0.
    y0 = 0.
    dx = round(18.e5 / nx)
    dy = dx
    h_mask = trues(nx,ny)
    grid = Grid(nx = nx, 
            ny = ny,   
            nσ = nσ, 
            x0 = x0, 
            y0 = y0, 
            dx = dx, 
            dy = dy,
            h_mask = h_mask)

    #bed
    function b(x,y)  
        B = 720 - 778.5 * x ./ (750e3)
        return B
    end

    #solver parameters
    maxiter_picard = 1
    solver_params = SolverParams(maxiter_picard = maxiter_picard)

    #Physical parameters
    default_thickness = 50.0 #set the initial condition this way
    accumulation_rate = 0.3
    sec_per_year = 365.25*24*60*60
    glen_a_ref = A * sec_per_year 

    params = Params(default_thickness = default_thickness, 
                    accumulation_rate = accumulation_rate,
                    glen_a_ref = glen_a_ref)
     
    #make the model
    model = Model(grid = grid,
    bed_elevation = b, 
    params = params, 
    solver_params = solver_params)

    
    #timestepping parameters
    niter0 = 0
    dt = 0.5
    end_time = 250.
    timestepping_params = TimesteppingParams(niter0 = niter0, 
                                            dt = dt, 
                                            end_time = end_time)

    folder = joinpath(@__DIR__, "./outputs_mismip_exp1")
    isdir(folder) && rm(folder, force = true, recursive = true)
    mkdir(folder) #make a clean folder for outputs
    outputs = (h = model.fields.gh.h,
                u = model.fields.gh.u,
                v = model.fields.gh.v,
                b = model.fields.gh.b,
                s = model.fields.gh.s,
                grfrac = model.fields.gh.grounded_fraction) #output velocities and thickness

    output_freq = 10.
    output_params = OutputParams(outputs = outputs, 
                            output_path = folder,
                            output_freq = output_freq,
                            output_format = "jld2",
                            zip_format = "none")

    simulation = Simulation(model = model, 
    timestepping_params = timestepping_params,
    output_params = output_params)

    #perform the simulation
    run_simulation!(simulation)

    return simulation
end

A = 4.6416e-24
simulation = MISMIP_EXP1_1(A)

################################################
function plot_evolution()
outfolder = joinpath(@__DIR__, "outputs_mismip_exp1")
files = [joinpath(outfolder, file) for file in readdir(outfolder) if endswith( joinpath(outfolder, file), ".jld2") ] 
nout = length(files)

#get the bed
dd = load(files[1])
bed  = dd["b"][:,1]
x    = dd["x"][:,1]
nx = length(bed[:,1])
h_out = zeros(nx, nout)
s_out = zeros(nx, nout)
ib_out = zeros(nx, nout)
t_out = zeros(1,nout)
 
#load solution files into matrix
for i = 1:nout
    fname = files[i]
    d = load(fname)
    
    h_out[:,i] = d["h"][:,1]
    ib_out[:,i] = d["s"][:,1] .- d["h"][:,1]
    s_out[:,i] = d["s"][:,1]
    t_out[i] = d["t"]
end
    
#plot things
plot1 = Plots.plot(x, bed, legend = false)
Plots.plot!(x, s_out)
Plots.plot!(x, ib_out)

display(plot1)
return nothing
end

plot_evolution()