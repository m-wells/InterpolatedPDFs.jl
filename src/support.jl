_get_typename(f::Function)::Core.TypeName = typeof(f).name
_get_mt(n::Core.TypeName)::Core.MethodTable = n.mt
_get_maxargs(mt::Core.MethodTable)::Int = mt.max_args
_get_maxargs(f::Function)::Int = _get_maxargs(_get_mt(_get_typename(f)))
_get_numargs(f::Function)::Int = _get_maxargs(f) - 1

function SupportFunction(f::Function, ::Type{T}, expr::Expr) where T<:Tuple
    fstr = string(prettify(expr))
    for t in fieldtypes(T)
        if isabstracttype(t)
            @warn """
                abstract types detected while constructing SupportFunction
                    f: $fstr
                    T: $T
                using abstract types carries a performance penality
                consider using concrete types
                """
            break
        end
    end
    n = length(fieldnames(T))
    N = _get_numargs(f)
    n == N || error("""
        Number of function arguments does not match number of argument types.
            f: $fstr
            _get_numargs(f) = $N
            typeof(args): $T
            length(fieldnames(T)) = $n
        """
       )
    return SupportFunction{T}(f, string(prettify(expr)))
end

function SupportFunction(f::Function, ::Type{T}, expr::Expr) where T
    N = nfields(methods(f))
    return SupportFunction(f, NTuple{N,T}, expr)
end

function SupportFunction(f::Function, expr::Expr)
    return SupportFunction(f, Any, expr)
end


"""
    @supportfunction expr T

Create a `SupportFunction` that evaluates an `expr` and takes arguments of type `T`.

# Examples
```julia-doctest
julia> sf = @supportfunction (x,y)->sin(x)*y>1 Tuple{Float64,Int}
SupportFunction:
   (x, y)->sin(x) * y > 1

julia> sf(π,3)
false

julia> sf(π/2,3)
true
```
"""
macro supportfunction(expr,T)
    esc(:(InterpolatedPDFs.SupportFunction($expr, $T, $(Expr(:quote, expr)))))
end

"""
    @supportfunction expr

Create a `SupportFunction` that evaluates an `expr` and takes arguments of type `Any`.
When possible use `@supportfunction expr T` where `T = Tuple{T1,T2,...}` where `T1,T2,...` are the types of the arguments.

# Examples
```julia-doctest
julia> sf = @supportfunction (x,y)->sin(x)*y>1 Tuple{Float64,Int}
SupportFunction:
   (x, y)->sin(x) * y > 1

julia> sf(π,3)
false

julia> sf(π/2,3)
true
```
"""
macro supportfunction(expr)
    esc(:(InterpolatedPDFs.SupportFunction($expr, $(Expr(:quote, expr)))))
end

#is_supported(sf::SupportFunction{T}, x::T) where T = sf(x...)
#function is_supported(sfs::Vector{SupportFunction{T}}, x::T) where T
#    return all(sf -> is_supported(sf,x), sfs)
#end

"""
    SupportedInterp(interp::Extrapolation, sfs::Vector{SupportFunction})

Create a supported interpolation.
"""
function SupportedInterp(interp::Extrapolation{Float64,N,ITP,IT,BC},
                         sfs::Vector{SupportFunction{S}}
                        ) where {N,ITP,IT,BC,S}
    knots = getknots(interp)
    # validate input
    for knot in Iterators.product(knots...)
        for sf in sfs
            # if knot is not supported and is not zero throw an error
            if !sf(knot) && !iszero(interp(knot...))
                throw(DomainError(knot, string(
                    "Found non-zero out of bound element.\n",
                    "knot = ", knot, "\n",
                    "sf(knot...) = ", sf(knot...), "\n",
                    "interp(knot...) = ", interp(knot...), "\n"
                   )))
            end
        end
    end

    binds = eachbin(coefficients(interp))
    boundary_mask = falses(size(binds))
    for (i,bind) in enumerate(binds)
        for sf in sfs
            if !all(sf, knots[bind]) && !all(x->!sf(x), knots[bind])
                boundary_mask[i] = true
            end
        end
    end

    return SupportedInterp{N,ITP,IT,BC,S}(interp, sfs, boundary_mask)
end

SupportedInterp(interp::Extrapolation, sf::SupportFunction) = SupportedInterp(interp, [sf])


"""
    SupportedInterp(knots::NTuple{N,AbstractVector}, coeffs::AbstractArray{T,N}, [sfs])

Convenience function that creates a `SupportedInterp`.
The optional `sfs` can be either a `SupportFunction` or a `Vector{SupportFunction}`.
"""
function SupportedInterp(knots::NTuple{N,AbstractVector},
                         coeffs::AbstractArray{T,N},
                         sfs = SupportFunction{Tuple{Nothing}}[]
                        ) where {N,T}
    interp = LinearInterpolation(knots, coeffs)
    return SupportedInterp(interp, sfs)
end

getknots(s::SupportedInterp) = getknots(s.interp)
coefficients(s::SupportedInterp) = coefficients(s.interp)

getknots(s::SupportedInterp, i) = getknots(s)[i]
coefficients(s::SupportedInterp, i) = coefficients(s.interp)[i]

isboundary(s::SupportedInterp, i) = s.boundary_mask[i]
