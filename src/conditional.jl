# here to catch cases where inds is already fully specified
sample(::Any, inds::NTuple{N,Int}) where N = error("""
    inds is fully specified unable to sample.
    inds = $(inds)
    """
   )

# here to catch cases where inds is already fully specified
samples(::Any, inds::NTuple{N,Int}, n::Int) where N = error("""
    inds is fully specified unable to sample.
    inds = $(inds)
    """
   )

"""
    indexrange(X::AbstractArray{T,N}, inds::Inds{N}, i::Int)

For dimension `i` of `X` return the range of possible indices specified by `inds`.
`inds` should be a tuple consisting of `Int`s and `Colon`s.

# Examples
```jldoctest
julia> using InterpolatedPDFs: indexrange

julia> indexrange(rand(4,3,2), (2,:,1), 2)
1:3

julia> indexrange(rand(4,3,2), (2,:,1), 3)
1:1

```
"""
function indexrange(x::AbstractArray{<:Any,N}, inds::T, i::Int) where {N,T<:Inds{N}}
    fieldtype(T,i) == Colon && return 1:size(x,i)
    ind = inds[i]::Int
    return ind:ind
end

"""
    sample(x::AbstractArray{Float64,N}, inds::Inds{N}) where N

Conditionally sample `x` using `inds`.

# Examples
```julia-repl
julia> x = rand(4,3,2); x./=sum(x);

julia> sample(x, (:,1,2))
CartesianIndex(2, 1, 2)

julia> sample(x, (3,:,2))
CartesianIndex(3, 3, 2)

```
"""
function sample(x::AbstractArray{Float64,N}, inds::Inds{N}) where N
    cinds = CartesianIndices(ntuple(i->indexrange(x, inds, i), Val(N)))
    GC.@preserve x begin
        let x = uview(x, inds...)
            vx = vec(x)
            sind = sample(vx)
            @inbounds return cinds[sind]
        end
    end
end

"""
    samples(x::AbstractArray{Float64,N}, inds::Inds{N}, n::Int) where N

Conditionally sample `x` `n` times using `inds`.

# Examples
```julia-repl
julia> x = rand(4,3,2); x./=sum(x);

julia> samples(x, (3,:,2), 2)
2-element Array{CartesianIndex{3},1}:
 CartesianIndex(3, 3, 2)
 CartesianIndex(3, 2, 2)

```
"""
function samples(x::AbstractArray{Float64,N}, inds::Inds{N}, n::Int) where N
    cinds = CartesianIndices(ntuple(i->indexrange(x, inds, i), Val(N)))
    GC.@preserve x begin
        let x = uview(x, inds...)
            vx = vec(x)
            sinds = samples(vx,n)
            @inbounds return cinds[sinds]
        end
    end
end
