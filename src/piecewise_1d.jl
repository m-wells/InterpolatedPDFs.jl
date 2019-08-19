"""
    ContinuousPiecewiseLinear{T,D} <: ContinuousUnivariateDistribution

A continuous univariate piecewise distribution.
The pdf, cdf, and inverse cdf are interpolated.
Using this construction directly requires the input to be normalized.

See also: [`fit_cpl`](@ref)

# Examples
```julia-repl
julia> x,y = [1.0, 2.0, 3.0], [0.75, 0.5, 0.25]
([1.0, 2.0, 3.0], [0.75, 0.5, 0.25])

julia> ContinuousPiecewiseLinear(x,y)
ContinuousPiecewiseLinear{Float64,1}(
pdf_itp: 3-element extrapolate(interpolate((::Array{Float64,1},), ::Array{Float64,1}, Gridded(Linear())), Throw()) with element type Float64:
 0.75
 0.5
 0.25
cdf_itp: 3-element extrapolate(interpolate((::Array{Float64,1},), ::Array{Float64,1}, Gridded(Linear())), Throw()) with element type Float64:
 0.0
 0.625
 1.0
invcdf_itp: 3-element extrapolate(interpolate((::Array{Float64,1},), ::Array{Float64,1}, Gridded(Linear())), Throw()) with element type Float64:
 1.0
 2.0
 3.0
)

julia> pdf(d,1.5)
0.625

julia> cdf(d,2.5)
0.8125

julia> quantile(d,0.5)
1.8
```
"""
struct ContinuousPiecewiseLinear{T,D} <: ContinuousUnivariateDistribution
    pdf_itp::Extrapolation
    cdf_itp::Extrapolation
    invcdf_itp::Extrapolation

    function ContinuousPiecewiseLinear(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Real}
        area = integrate(x,y)
        area ≈ 1 || error("Input is not normalized. integrate(x,y) = ", area)
        issorted(x) || error("x is not in ascending order")

        pdf_itp = LinearInterpolation(x,y)

        cdf_y = cumul_integrate(x,y)
        cdf_itp = LinearInterpolation(x,cdf_y)
        invcdf_itp = LinearInterpolation(cdf_y,x)

        new{T,1}(pdf_itp,cdf_itp,invcdf_itp)
    end
end

pdf(d::ContinuousPiecewiseLinear, x::T) where {T<:Real} = d.pdf_itp(x)
cdf(d::ContinuousPiecewiseLinear, x::T) where {T<:Real} = d.cdf_itp(x)
quantile(d::ContinuousPiecewiseLinear, x::T) where {T<:Real} = d.invcdf_itp(x)

pdf(d::ContinuousPiecewiseLinear, x::AbstractVector{T}) where {T<:Real} = d.pdf_itp(x)
cdf(d::ContinuousPiecewiseLinear, x::AbstractVector{T}) where {T<:Real} = d.cdf_itp(x)
quantile(d::ContinuousPiecewiseLinear, x::AbstractVector{T}) where {T<:Real} = d.invcdf_itp(x)

midpoints(x::AbstractRange) = range(first(x)+step(x)*(1//2), last(x)-step(x)*(1//2), length=length(x)-1)
midpoints(x::AbstractVector) = [(x[i]+x[i+1])*(1//2) for i in eachindex(x)[1:end-1]]

"""
    fit_cpl(x::AbstractArray, s::AbstractArray)

Fit s with a ContinuousPiecewiseLinear distribution using x for the breakpoints.

See also: [`ContinuousPiecewiseLinear`](@ref)

# Examples
```julia-repl
julia> x = range(0,pi,length=5);

julia> s = acos.(rand(100));

julia> d = fit_cpl(x,s)
ContinuousPiecewiseLinear{Float64,1}(
pdf_itp: 5-element extrapolate(scale(interpolate(::Array{Float64,1}, BSpline(Linear())), (0.0:0.39269908169872414:1.5707963267948966,)), Throw()) with element type Float64:
 0.15242714502015237
 0.36726001464827923
 0.7290336549033926
 0.8405744820809898
 1.0667947306551753
cdf_itp: 5-element extrapolate(scale(interpolate(::Array{Float64,1}, BSpline(Linear())), (0.0:0.39269908169872414:1.5707963267948966,)), Throw()) with element type Float64:
 0.0
 0.10204033518620566
 0.31729709383873866
 0.6254889308490448
 1.0
invcdf_itp: 5-element extrapolate(interpolate((::Array{Float64,1},), ::Array{Float64,1}, Gridded(Linear())), Throw()) with element type Float64:
 0.0
 0.39269908169872414
 0.7853981633974483
 1.1780972450961724
 1.5707963267948966
)
"""
function fit_cpl(x::AbstractArray, s::AbstractArray)
    !issorted(x) && error("breakpoints are not in ascending order")
    minimum(s) < x[1] && error("encounter values below minimum breakpoint")
    maximum(s) > x[end] && error("encounter values above maximum breakpoint")

    length(x) < 2 && error("need a minimum of 2 breakpoints")

    cdf_x = [first(x), sort(s)..., last(x)]
    cdf_y = range(0,1,length=length(cdf_x))
    cdf_itp = LinearInterpolation(cdf_x,cdf_y)

    invcdf_itp = LinearInterpolation(cdf_y,cdf_x)

    xmid = [first(x), midpoints(x)..., last(x)]
    #Δxmid = diff(xmid)
    ymid = cdf_itp(xmid)
    #Δymid = diff(ymid)

    #p = similar(x)
    #p[1] = first(ymid)/(first(xmid) - first(x))
    #p[end] = (1 - last(ymid))/(last(x) - last(xmid))
    #p[2:end-1] = Δymid./Δxmid
    p = diff(ymid)./diff(xmid)

    return ContinuousPiecewiseLinear(x, p)
end

using Random
import Base.rand
rand(rng::AbstractRNG, d::ContinuousPiecewiseLinear) = d.invcdf_itp(rand(rng))
rand(rng::AbstractRNG, d::ContinuousPiecewiseLinear, n::Int) = d.invcdf_itp(rand(rng,n))
