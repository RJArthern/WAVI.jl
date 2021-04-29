function driver(fname,end_time, output_time, exist_flag)

    # MISMIP+
    nx = 320
    ny = 40
    nσ = 4
    x0 = 0.0
    y0 = -40000.0
    dx = 2000.0
    dy = 2000.0
    dt = 0.025
    xx=[x0+(i-0.5)*dx for i=1:nx, j=1:ny]
    yy=[y0+(j-0.5)*dy for i=1:nx, j=1:ny]
    accumulation_rate=0.3
    
    if exist_flag
    h = ncread(fname, "h");
    h = h[:,:,end]
    starting_thickness=h
    print("Found an input file, picking up /n")
    else
    starting_thickness=100.0.*ones(nx,ny)
    print("Did not find an input file, starting afresh /n")
    end
    
    #Model domain mask on h-grid
    h_mask=trues(nx,ny)
    
    #Homogenous Dirichlet boundary conditions
    u_iszero=falses(nx+1,ny)
    u_iszero[1,:].=true
    v_iszero=falses(nx,ny+1)
    v_iszero[:,1].=true
    v_iszero[:,end].=true
    
    #MISMIP+ bed elevation
    function bed_elev_function(x,y)
    xbar = 300000.0
    b0 = -150.0; b2 = -728.8; b4 = 343.91; b6 = -50.75
    wc = 24000.0; fc = 4000.0; dc = 500.0
    bx(x)=b0+b2*(x/xbar)^2+b4*(x/xbar)^4+b6*(x/xbar)^6
    by(y)= dc*( (1+exp(-2(y-wc)/fc))^(-1) + (1+exp(2(y+wc)/fc))^(-1) )
    b = max(bx(x) + by(y), -720.0)
    return b
    end
    bed_elevation=bed_elev_function.(xx,yy)
    
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
    print("Starting WAVI run...")
    wavi=start(params)
    run!(wavi)
    
    ########### intialize the netcdf file ###########
    if exist_flag
    tcurr = ncread(fname, "t");
    tcurr = tcurr[end]
    else
    tcurr = 0.
    init_nc(fname,wavi)
    end
    
    ########## Do the timestepping #########
    while tcurr < end_time
        run!(wavi)
        tcurr = tcurr + params.dt
    #    print("\r")
    #    print("Completed: ",round(tcurr), " years")
       
       if abs(mod(tcurr, output_time)) < params.dt 
            update_nc!(fname, wavi, tcurr)
        print("updating nc file...t = ")
        print(tcurr)
       end
    end
    print("Finished simulation")
    return wavi
    end
    
    ###########################################################
    using WAVI, WAVItools, NetCDF
    fname = "MISMIP_013.nc"
    end_time = 10000 #endtime of the simulation in yrs
    output_time = 100 #
    
    #work out whether filename exists or not:
    file_exists = isfile(fname) #returns true if fname exists. In which case, run from this.
    
    driver(fname, end_time, output_time, file_exists)
    