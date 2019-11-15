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

function Base.getindex(knots::NTuple{N,AbstractVector}, cind::CartesianIndex{N}) where {T,N}
    all(i->(1 ≤ i ≤ length(knots[i])), 1:N) || throw(BoundsError(knots, cind))
    return getindex.(knots,cind.I)
end

function Base.getindex(knots::NTuple{N,AbstractVector}, cinds::CartUnitInds{N}) where N
    return map(x->getindex(knots,x), cinds)
end


"""
    randpnt(x1::Real, x2::Real)

Return a uniformly random value between `x1` and `x2`.
"""
randpnt(x1::Real,x2::Real) = x1 + rand()*(x2-x1)

"""
    randpnt(x::AbstractVector, inds::UnitRange{Int})

Return a random value between `x[first(inds)]` and `x[last(inds)]`.
"""
function randpnt(x::AbstractVector, inds::UnitRange{Int})
    i1,i2 = extrema(inds)
    @inbounds return randpnt(x[i1],x[i2])
end

"""
    randpnt(knots::NTuple{N,AbstractVector}, cinds::CartUnitInds{N})

Generate a uniformly random point (`Tuple`) for a given bin.
The corners of the bin are given by `cinds` and `knots` define the coordinate axes.

# Examples
```julia-repl
julia> knots = ([0.0,0.2,0.4,0.6],[1.0,2.0,3.0],[0,π]);

julia> inds = CartesianIndices((3:4,2:3,1:2));

julia> randpnt(knots,inds)
(0.5376059603415213, 2.2621926859178796, 0.2515135019823618)

```
"""
function randpnt(knots::NTuple{N,AbstractVector}, cinds::CartUnitInds{N}) where N
    @inbounds return ntuple(i -> randpnt(knots[i], cinds.indices[i]), Val(N))
end

function randpnt(s::SupportedInterp{N,ITP,IT,BC,S},
                 cinds::CartUnitInds{N}
                ) where {N,ITP,IT,BC,S}
    knots = getknots(s)
    return randpnt(knots, cinds)
end
