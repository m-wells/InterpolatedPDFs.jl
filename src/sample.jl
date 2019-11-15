"""
    sample(x::AbstractVector, [sumx])

Returns a randomly selected index according to the values of `x`.
The elements of `x` are only assumed to be positive and are not required to sum to one.
A precomputed sum can be given via `sumx` otherwise one with be computed.
This does not perform any allocation.
"""
function sample(x::AbstractVector{<:Real}, sumx::Real = sum(x))
    c = 0.0
    v = sumx*rand()
    for (i,p) in enumerate(x)
        c+=p
        c ≥ v && return i
    end
    error("unable to draw index, check input => ", x)
end

"""
    sample(x::AbstractArray, [sumx])

Returns a randomly selected CartesianIndex according to the values of `x`.
The elements of `x` are only assumed to be positive and are not required to sum to one.
A precomputed sum can be given via `sumx` otherwise one with be computed.
This does not perform any allocation.
"""
function sample(x::AbstractArray{<:Real}, sumx::Real = sum(x))
    cinds = CartesianIndices(x)
    @uviews x begin
        vx = vec(x)
        ind = sample(vx, sumx)
        @inbounds cinds[ind]
    end
end

"""
    samples(x::AbstractVector, n::Int)

Returns `n` randomly selected indices according to the values of `x`.
The elements of `x` are only assumed to be positive and are not required to sum to one.
This returns a `Vector{Int}`
"""
function samples(x::AbstractVector{<:Real}, n::Int)
    c = cumsum(x)
    v = last(c)*rand(n)

    retval = Vector{Int}(undef,n)
    for i in eachindex(retval)
        @inbounds retval[i] = searchsortedfirst(c,v[i])
    end

    return retval
end

"""
    samples(x::AbstractArray{T,N}, n::Int) where {T,N}

Returns `n` randomly selected CartesianIndices according to the values of `x`.
The elements of `x` are only assumed to be positive and are not required to sum to one.
This returns `Vector{CartesianIndex{N}}`
"""
function samples(x::AbstractArray{T,N}, n::Int) where {T<:Real,N}
    @uviews x begin
        vx = vec(x)
        sample_inds = samples(vx,n)
        cinds = CartesianIndices(x)
        @inbounds cinds[sample_inds]
    end
end

"""
    sample(p::PMF)

Returns a randomly selected bin from the supplied probability mass function (`PMF`).
The bin is defined by the CartesianIndices of the corners.
"""
function sample(pmf::PMF{N,T}) where {N,T}
    cind = sample(pmf.vals, 1.0)
    @inbounds return CartesianIndices(ntuple(i -> size(pmf,i) > 1 ?  (cind[i]:cind[i]+1) : (1:1),
                                             Val(N)
                                            )
                                     )
end

#"""
#    samples(p::PMF, n::Int)
#
#Returns a randomly selected index from the supplied probability mass function (`PMF`).
#"""
#samples(pmf::PMF, n::Int) = samples(pmf.vals, n)
#
#"""
#    sample(p::PMF, inds::Inds)
#
#Conditionally sample the probability mass function `p` using `inds`.
#"""
#sample(pmf::PMF{N,<:AbstractArray}, inds::Inds{N}) where N = sample(pmf.vals, inds)
#
#"""
#    samples(p::PMF, inds::Inds, n::Int) where {T,N}
#
#Conditionally sample the probability mass function `p` `n` times using `inds`.
#"""
#samples(pmf::PMF{N,<:AbstractArray}, inds::Inds{N}, n::Int) where N = sample(pmf.vals, inds, n)

sample(s::SupportedInterp) = error("sampling a SupportedInterp requires specifying a bin")

# this performs much better than maximum(x[inds])
function _max(x,inds)
    xmax = -Inf
    for i in inds
        x[i] > xmax && (xmax = x[i])
    end
    return xmax
end

"""
    sample(s::SupportedInterp, bin::CartUnitInds)

Sample a specified bin from the `SupportedInterp`.
"""
function sample(s::SupportedInterp{N,ITP,IT,BC,S},
                bin::CartUnitInds{N}
               ) where {N,ITP,IT,BC,S}
    maxp = _max(s, bin)
    while maxp > 0
        pnt = randpnt(s, bin)
        s(pnt) < rand()*maxp && return pnt
    end
    throw(DomainError(bin,"""
                      No suitable coefficients found!
                      s[bin] = $(s[bin])
                      """
                     )
         )
end

"""
    sample(m::MultivariatePDF)

Sample from a multivariate pdf.
"""
function sample(m::MultivariatePDF)
    bin = sample(m.pmf)
    return sample(m.pdf, bin)
end
