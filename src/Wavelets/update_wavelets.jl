export update_wavelets!, wavelet_matrix

"""
    update_wavelets!(model::AbstractModel)

Compute wavelet transform of velocities to define the coarse grid used in multigrid preconditioner.
"""
function update_wavelets!(model::AbstractModel{T}) where {T}
    @unpack wu,wv,gu,gv=model.fields
    @unpack params,solver_params=model
    s = stencil_scratch!(model)

    copyto!(wu.wavelets, gu.u)
    @inbounds for i in eachindex(wu.wavelets, gu.mask)
        wu.wavelets[i] *= gu.mask[i]
    end
    haar_dwt!(wu.wavelets, s.haar_u_tmp, wu.levels)

    copyto!(wv.wavelets, gv.v)
    @inbounds for i in eachindex(wv.wavelets, gv.mask)
        wv.wavelets[i] *= gv.mask[i]
    end
    haar_dwt!(wv.wavelets, s.haar_v_tmp, wv.levels)

    wu.mask .= (abs.(wu.wavelets) .>= solver_params.wavelet_threshold)
    wv.mask .= (abs.(wv.wavelets) .>= solver_params.wavelet_threshold)

    wu.n[] = count(wu.mask)
    wv.n[] = count(wv.mask)

    wu.correction_coarse[] = zeros(T,wu.n[])
    wv.correction_coarse[] = zeros(T,wv.n[])

    fill_index_map!(wu.index_map, wu.mask)
    fill_index_map!(wv.index_map, wv.mask)

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