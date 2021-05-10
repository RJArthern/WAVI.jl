"""
    update_wavelets!(model::AbstractModel)

Compute wavelet transform of velocities to define the coarse grid used in multigrid preconditioner.
"""
function update_wavelets!(model::AbstractModel)
    @unpack wu,wv,gu,gv,params,solver_params=model

    wu.wavelets[:] .= gu.dwt*(gu.crop*gu.u[:])
    wv.wavelets[:] .= gv.dwt*(gv.crop*gv.v[:])

    wu.mask .= (abs.(wu.wavelets) .>= solver_params.wavelet_threshold)
    wv.mask .= (abs.(wv.wavelets) .>= solver_params.wavelet_threshold)

    wu.n[] = count(wu.mask)
    wv.n[] = count(wv.mask)

    wu.crop[] .= Diagonal(float(wu.mask[:]))
    wv.crop[] .= Diagonal(float(wv.mask[:]))

    wu.samp[] = wu.crop[][wu.mask[:],:]
    wv.samp[] = wv.crop[][wv.mask[:],:]

    wu.spread[] = sparse(wu.samp[]')
    wv.spread[] = sparse(wv.samp[]')

    return model
end
"""
    wavelet_matrix(n,levels,direction)

Compute matrix used to apply wavelet transform in one dimension.
"""
function wavelet_matrix(n,levels,direction);

    identity=spdiagm(n,n, 0 => ones(n))

    isforward = lowercase(direction) == "forward"
    isreverse = lowercase(direction) == "reverse"

    if isforward
        fwav=identity
    elseif isreverse
        rwav=identity
    else
        error("Direction of wavelet transform not defined properly")
    end

    step=1

    for j=levels:-1:0

        step=step*2

        ixodd=1:step:(n-div(step,2))
        ixeven=ixodd.+div(step,2)
        ixothers=setdiff(1:n, union(ixodd,ixeven))

        odds=identity[ixodd,:]
        evens=identity[ixeven,:]
        others=identity[ixothers,:]

        if isforward
            #Forward Haar wavelet lifting scheme
            fwav = (others'*others+0.5*(evens'*(evens - odds) + odds'*(odds + evens)))*fwav
        elseif isreverse
            #Reverse Haar wavelet lifting scheme
            rwav = rwav*(others'*others + evens'*(evens + odds) + odds'*(odds - evens))
        else
            error("Direction of wavelet transform not defined properly")
        end

    end

    if isforward
            #Forward Haar wavelet lifting scheme
            return fwav
    elseif isreverse
            #Reverse Haar wavelet lifting scheme
            return rwav
    else
            error("Direction of wavelet transform not defined properly")
    end

end