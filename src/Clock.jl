#mutable clock structure to store time info
mutable struct Clock{T <: Real, N <: Integer}
    n_iter::N
    time::T
end

#clock constructor
function Clock(;
                n_iter = 0,
                time = 0)
    return Clock(n_iter, time)
end