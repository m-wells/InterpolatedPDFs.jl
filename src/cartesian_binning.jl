"""
    EachBin(dims::Dims)
    EachBin(X::AbstractArray)

An iterable `EachBin` object that returns the CartesianIndices of the corners for each bin of the given `dims` or for an `AbstractArray`.
"""
struct EachBin{N}
    indices::CartesianIndices{N,NTuple{N,UnitRange{Int64}}}
    corners::CartesianIndices{N,NTuple{N,UnitRange{Int64}}}
end

function Base.show(io::IO, b::EachBin)
    println(io, "EachBin.indices => ", b.indices)
    println(io, "EachBin.corners => ", b.corners)
end
function Base.show(io::IO, ::MIME"text/plain", b::EachBin)
    print(io, "EachBin.indices => ")
    show(io, MIME("text/plain"), b.indices)
    print(io, "\nEachBin.corners => ")
    show(io, MIME("text/plain"), b.corners)
end

"""
    eachbin(X::AbstractArray)

Create an iterable `EachBin` object that yields the CartesianIndices of the corners of each bin in an array `X`.
If knots are a MxN array then there will be (M-1)x(N-1) bins.
Dimensions with a single knot are allowed and supported.
"""
function eachbin(dims::Dims{N}) where N
    csz = ntuple(i -> dims[i] > 1 ? dims[i] - 1 : dims[i], Val(N))
    corners = CartesianIndices(ntuple(i -> csz[i] > 1 ? (0:1) : (0:0), Val(N)))
    return EachBin{N}(CartesianIndices(csz), corners)
end
eachbin(arr::AbstractArray) = eachbin(size(arr))
eachbin(s::SupportedInterp) = eachbin(coefficients(s))

Base.length(b::EachBin) = length(b.indices)
Base.size(b::EachBin) = size(b.indices)

function Base.iterate(b::EachBin, state=1)
    state > length(b) && return nothing
    return (b.indices[state] .+ b.corners, state+1)
end

function Base.getindex(b::EachBin, i::Int)
    1 ≤ i ≤ length(b) || throw(BoundsError(b, i))
    return b.indices[i] .+ b.corners
end

Base.firstindex(b::EachBin) = firstindex(b.indices)
Base.lastindex(b::EachBin) = lastindex(b.indices)

function Base.getindex(knots::NTuple{N,AbstractVector}, cind::CartesianIndex{N}) where {T,N}
    all(i->(1 ≤ i ≤ length(knots[i])), 1:N) || throw(BoundsError(knots, cind))
    return getindex.(knots,cind.I)
end

function Base.getindex(knots::NTuple{N,AbstractVector},
                       cinds::CartesianIndices{N,NTuple{N,UnitRange{Int}}}
                      ) where N
    return map(x->getindex(knots,x), cinds)
end
