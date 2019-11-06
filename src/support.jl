_get_typename(f::Function)::Core.TypeName = typeof(f).name
_get_mt(n::Core.TypeName)::Core.MethodTable = n.mt
_get_maxargs(mt::Core.MethodTable)::Int = mt.max_args

_get_maxargs(f::Function)::Int = _get_maxargs(_get_mt(_get_typename(f)))
_get_numargs(f::Function)::Int = _get_maxargs(f) - 1

"""
    SupportFunction(f::Function, expr::String)

Stores function and a human readable version of the function.
"""
struct SupportFunction{T}
    f::FunctionWrapper{Bool,T}
    expr::String

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
        return new{T}(f, string(prettify(expr)))
    end

    function SupportFunction(f::Function, ::Type{T}, expr::Expr) where T
        N = nfields(methods(f))
        return SupportFunction(f, NTuple{N,T}, expr)
    end

    function SupportFunction(f::Function, expr::Expr)
        return SupportFunction(f, Any, expr)
    end
end

"""
    @supportfunction expr T

Create a `SupportFunction` that evaluates an `expr` and takes arguments of type `T`.

# Examples
```julia-repl
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

macro supportfunction(expr)
    esc(:(InterpolatedPDFs.SupportFunction($expr, $(Expr(:quote, expr)))))
end

(bf::SupportFunction)(args...) = bf.f(args...)::Bool
(bf::SupportFunction)(args::Tuple) = (bf)(args...)

Base.show(io::IO, sf::SupportFunction) = print(io, sf.expr)
Base.show(io::IO, ::MIME"text/plain", sf::SupportFunction) =
           print(io, "SupportFunction:\n   ", sf.expr)

"""
Wrapper around an interpolation that allows for well defined boundary conditions
"""
struct SupportedInterp{T,N,ITP,IT,BC,S}
    interp::Extrapolation{T,N,ITP,IT,BC}
    support_functions::Vector{SupportFunction{S}}
    boundary_mask::BitArray{N}

    function SupportedInterp(interp::Extrapolation{T,N,ITP,IT,BC},
                             sfs::Vector{SupportFunction{S}} = SupportFunction{Tuple{Nothing}}[]
                            ) where {T,N,ITP,IT,BC,S}
        knots = getknots(interp)
        # validate input
        for knot in Iterators.product(knots...)
            for sf in sfs
                # if knot is not supported and is not zero throw an error
                if !sf(knot...) && !iszero(interp(knot...))
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

        return new{T,N,ITP,IT,BC,S}(interp, sfs, boundary_mask)
    end

    function SupportedInterp(interp::Extrapolation, sf::SupportFunction)
        return SupportedInterp(interp, [sf])
    end

    function SupportedInterp(knots, coeffs,
                             sfs::Vector{SupportFunction{T}} = SupportFunction{Tuple{Nothing}}[]
                            ) where T
        interp = LinearInterpolation(knots, coeffs)
        return SupportedInterp(interp, sfs)
    end
    
    function SupportedInterp(knots, coeffs, sf::SupportFunction)
        interp = LinearInterpolation(knots, coeffs)
        return SupportedInterp(interp, [sf])
    end
end

function (sinterp::SupportedInterp)(args...)
    for f in sinterp.support_functions
        # if any of these evaluate false then return 0.0
        !f(args...) && return 0.0
    end
    # otherwise use the interpolation
    return sinterp.interp(args...)
end
(s::SupportedInterp)(args::Tuple) = (s)(args...)

function Base.getindex(sinterp::SupportedInterp{T,N,ITP,IT,BC,S},
                       cinds::CartesianIndices{N,NTuple{N,UnitRange{Int}}}
                      ) where {T,N,ITP,IT,BC,S}
    return coefficients(sinterp)[cinds]
end

function Base.show(io::IO, s::SupportedInterp)
    println(io, s.interp)
    println(io, s.support_functions)
    print(io, s.boundary_mask)
end

function Base.show(io::IO, ::MIME"text/plain", s::SupportedInterp)
    println(io, "SupportedInterp")
    print(io, "interp => ")
    show(io, MIME("text/plain"), s.interp)
    print(io, "\nsupport_functions => ")
    show(io, MIME("text/plain"), s.support_functions)
    print(io, "\nboundary_mask => ")
    show(io, MIME("text/plain"), s.boundary_mask)
end

getknots(s::SupportedInterp) = getknots(s.interp)
coefficients(s::SupportedInterp) = coefficients(s.interp)

getknots(s::SupportedInterp, cind) = getknots(getknots(s), cind)
coefficients(s::SupportedInterp, cind) = coefficients(s.interp)[cind]

isboundary(s::SupportedInterp, cind::CartesianIndex) = s.boundary_mask[cind]
isboundary(s::SupportedInterp, ind::Int) = s.boundary_mask[ind]
