struct BinfileMeltRate{PC, M, N <: Int} <: AbstractMeltRateModel{PC,M}
    melt_partial_cell::PC     #specify whether melt applies to partial cells or not
    input_filename::String    #specify melt filename
    nx::N                     #number of grid points in x 
    ny::N                     #number of grid points in y
    niter_read::N             #specify how often to read the input file
    melt_rate::M              #stores the melt rate
end

function BinfileMeltRate(;
                        input_filename = nothing,
                        melt_partial_cell = true,
                        niter_read = true,
                        nx = nothing,
                        ny = nothing)

    #input file exceptions
    ~(input_filename === nothing) || throw(ArgumentError("You must pass an input filename"))
    isfile(input_filename) || throw(ArgumentError("Did not find the input file specified in InputFileMeltRateModel"))
    ((nx === nothing) || (ny === nothing)) || throw(ArgumentError("You must pass grid co-ordinates nx and ny to BinfileMeltRate"))
    
    

    return InputFileMeltRateModel(melt_partial_cell, input_filename, niter_read, melt_rate)
end


    
