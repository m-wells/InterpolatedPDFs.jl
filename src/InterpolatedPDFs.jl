module InterpolatedPDFs

#using Distributions
#using NumericalIntegration

using FunctionWrappers: FunctionWrapper
using MacroTools
using Base.Iterators

using Interpolations; using Interpolations: Extrapolation
import Interpolations: coefficients, getknots
coefficients(etp::Extrapolation) = coefficients(etp.itp)

const IntOrColon = Union{Int,Colon}
const Inds{N} = Tuple{Vararg{<:IntOrColon,N}}
const CartUnitInds{N} = CartesianIndices{N,NTuple{N,UnitRange{Int}}}

using UnsafeArrays

export sample, samples, @supportfunction, SupportedInterp, MultivariatePDF
#export LinearInterpolatedPDF, pdf, cdf, quantile

#include("./linear_1d.jl")

"""
    EachBin(dims::Dims)
    EachBin(X::AbstractArray)

An iterable `EachBin` object that returns the CartesianIndices of the corners for each bin of the given `dims` or for an `AbstractArray`.
"""
struct EachBin{N}
    indices::CartUnitInds{N}
    corners::CartUnitInds{N}
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

############################################################################################
############################################################################################
############################################################################################

"""
    SupportFunction(f::Function, expr::String)

Stores function and a human readable version of the function.
"""
struct SupportFunction{T}
    f::FunctionWrapper{Bool,T}
    expr::String
end

Base.show(io::IO, sf::SupportFunction) = print(io, sf.expr)
function Base.show(io::IO, ::MIME"text/plain", sf::SupportFunction)
    print(io, "SupportFunction:\n   ", sf.expr)
end

(sf::SupportFunction{S})(x::Vararg{T,N}) where {T,N,S<:NTuple{N,T}} = sf.f(x...)
(sf::SupportFunction{T})(x::T) where T = sf.f(x...)

############################################################################################
############################################################################################
############################################################################################

"""
    SupportedInterp{N,ITP,IT,BC,S}

Wrapper around an interpolation that allows for well defined boundary conditions
"""
struct SupportedInterp{N,ITP,IT,BC,S} <:AbstractArray{Float64,N}
    interp::Extrapolation{Float64,N,ITP,IT,BC}
    support_functions::Vector{SupportFunction{S}}
    boundary_mask::BitArray{N}
end

function Base.show(io::IO, s::SupportedInterp)
    show(io, s.interp)
    show(io, s.support_functions)
    show(io, s.boundary_mask)
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

"""
    (sinterp::SupportedInterp)(args...)

`SupportedInterp` is a callable object that will perform interpolation within the provided support.
"""
function (sinterp::SupportedInterp{N,ITP,IT,BC,S})(args::Vararg{Real,N}) where {N,ITP,IT,BC,S}
    for f in sinterp.support_functions
        # if any of these evaluate false then return 0.0
        f(args...) || return 0.0
    end
    # otherwise use the interpolation
    return sinterp.interp(args...)
end
(s::SupportedInterp)(args::Tuple) = (s)(args...)

Base.size(s::SupportedInterp) = size(s.interp)
Base.getindex(s::SupportedInterp, i) = getindex(coefficients(s), i)

############################################################################################
############################################################################################
############################################################################################

"""
    PMF(x::AbstractArray)

Construct a Probability Mass Function.
The input array `x` must contain only non-negative elements and must sum to 1.
"""
struct PMF{N,T<:AbstractArray{Float64,N}} <: AbstractArray{Float64,N}
    vals::T

    function PMF(vals::T) where {N,T<:AbstractArray{Float64,N}}
        isapprox(sum(vals), 1) || error("""
            Input does not sum to 1!
            sum(input) = $(sum(vals))
            """
           )
        any(signbit, vals) && error("""
            Input should only contain non-negative elements!
            """
           )
        new{N,T}(vals)
    end
end

Base.show(io::IO, pmf::PMF) = show(io::IO, pmf.vals)
function Base.show(io::IO, ::MIME"text/plain", pmf::PMF)
    show(io::IO, MIME("text/plain"), pmf.vals)
end

Base.size(p::PMF) = size(p.vals)
Base.getindex(p::PMF, i) = getindex(p.vals, i)

############################################################################################
############################################################################################
############################################################################################

"""
    MultivariatePDF

A grid approximated multivariate pdf.
"""
struct MultivariatePDF{N,ITP,IT,BC,S}
    pdf::SupportedInterp{N,ITP,IT,BC,S}
    pmf::PMF{N,Array{Float64,N}}
end

############################################################################################
############################################################################################
############################################################################################

include("./cartesian_binning.jl")
include("./support.jl")
include("./multivariate.jl")
include("./sample.jl")
include("./conditional.jl")

end # module
