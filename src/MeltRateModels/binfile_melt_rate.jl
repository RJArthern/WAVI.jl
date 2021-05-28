struct BinfileMeltRate{PC, M, N <: Int} <: AbstractMeltRateModel{PC,M}
    melt_partial_cell::PC     #specify whether melt applies to partial cells or not
    input_filename::String    #specify melt filename
    niter_read::N             #specify how often to read the input file
    melt_rate::M              #stores the melt rate
end

function BinfileMeltRate(;
                        input_filename = nothing,
                        melt_partial_cell = true,
                        niter_read = 1,
                        nx = nothing,
                        ny = nothing)

    #input file exceptions
    ~(input_filename === nothing) || throw(ArgumentError("You must pass an input filename"))
    isfile(input_filename) || throw(ArgumentError("Did not find the specified binary file"))
    (~(nx === nothing) && ~(ny === nothing)) || throw(ArgumentError("You must pass grid co-ordinates nx and ny to BinfileMeltRate"))
    
    #check input file size compatibility
    file_size = stat(input_filename).size #bytes in input file must match matrix returned
    empty_matrix = zeros(nx,ny)
    try 
        read!(input_filename, empty_matrix)
    catch
        Error("Input file read error")
    end
    (file_size == sizeof(empty_matrix)) || throw(DimensionMismatch("Size of input file incompatible with specified nx, ny"))

    #assemble the matrix
    melt_rate = zeros(nx,ny)

    return BinfileMeltRate(melt_partial_cell, input_filename, niter_read, melt_rate)
end


    

function update_melt_rate_model!(melt_model::BinfileMeltRate, model)
    @unpack melt_rate, input_filename = melt_model
    #read the input file
    read!(input_filename, melt_rate)
    melt_rate .= hton.(melt_rate)

#    if melt_model.melt_partial_cell
#        melt_rate .= melt_rate.* (1 .- model.fields.gh.grounded_fraction)
#    else
#        grounded_matrix =  model.fields.gh.grounded_fraction .- mod.(model.fields.gh.grounded_fraction,1) #send any partially floating cells to zero grounded fraction
#        melt_rate .= melt_rate.* (1 .- grounded_matrix)
#    end
    return nothing
end