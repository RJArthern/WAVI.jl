struct BinfileMeltRate <: AbstractMeltRate
    input_filename::String    #specify melt filename
end

function BinfileMeltRate(;
                        input_filename = nothing)

    #input file exceptions
    ~(input_filename === nothing) || throw(ArgumentError("You must pass an input filename"))
    isfile(input_filename) || throw(ArgumentError("Did not find the specified binary file"))
    
    return BinfileMeltRate(input_filename)
end


function update_melt_rate!(melt_model::BinfileMeltRate, fields, grid, clock)
    @unpack input_filename = melt_model
    @unpack basal_melt = fields.gh
    
    file_size = stat(input_filename).size #bytes in input file must match matrix returned
    melt_rate = zeros(grid.nx,grid.ny)
    (file_size == sizeof(melt_rate)) || throw(DimensionMismatch("Size of input file incompatible with specified nx, ny"))
    try 
        read!(input_filename, melt_rate)
    catch
        Error("Input file read error")
    end

   # melt_rate .= ntoh.(melt_rate)
    basal_melt[:] = melt_rate[:]
    return nothing
end

