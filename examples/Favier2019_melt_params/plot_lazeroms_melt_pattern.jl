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
    pme = PlumeEmulator()
    add_melt_rate_model!(model,pme) #this updates the melt rate in model

    return model
end



model = plot_lazeroms_melt_pattern()


















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