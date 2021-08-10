using WAVI 
#using PyPlot
#using LinearAlgebra
function plot_lazeroms_melt_pattern()
    #Grid and boundary conditions
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

    #solver parameters
    maxiter_picard = 1
    solver_params = SolverParams(maxiter_picard = maxiter_picard)

    #Physical parameters
    default_thickness = 100.0 #set the initial condition this way
    accumulation_rate = 0.3
    params = Params(default_thickness = default_thickness, 
                    accumulation_rate = accumulation_rate)
    

    #Inputing thickness profile, reading from binary file 
    fname = joinpath(dirname(@__FILE__), "data",  "MISMIP_ice0_2km_SteadyThickness.bin")
    h = Array{Float64, 2}(undef, nx, ny)
    read!(fname,h)
    h = ntoh.(h)        


    #make the model
    initial_conditions = InitialConditions(initial_thickness = h) #Defining initial condition of ice sheet
    model = Model(grid = grid,
                     bed_elevation = bed, 
                     params = params, 
                     solver_params = solver_params,
                     initial_conditions = initial_conditions)


    #embed the ice model with melt rate model
    #melt_rate_model = QuadraticMeltModel(model.fields.gh.h,
    #melt_partial_cell = true)
    #add_melt_rate_model!(model,melt_rate_model) #this updates the melt rate in model

    return h
end

function quad_melt_plot()
#@time simulation = MISMIP_PLUS();
sim = MISMIP_PLUS()
x = range(0, stop= 40*2., length = 40) #In km
y = range(-40., stop = 320*2., length = 320)

xgrid = repeat(x', 320 ,1)
ygrid = repeat(y, 1, 40)

basal_melt = sim.model.fields.gh.basal_melt


fig = figure()
cp = contourf(xgrid, ygrid, basal_melt, linewidth=2.0)
contour(xgrid,ygrid,sim.model.fields.gh.grounded_fraction)
PyPlot.ylim(400., 320*2.)
colorbar(cp)
PyPlot.title("Contour Plot of ice base melt rate m/yr")
savefig("plot.png")
display(fig)
return nothing
end

plot_lazeroms_melt_pattern()