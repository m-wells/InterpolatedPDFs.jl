"""
    LinearInterpolatedPDF{T,N,ITP,IT} <: ContinuousUnivariateDistribution

A continuous univariate linearly interpolated distribution.
The pdf, cdf, and inverse cdf are interpolated.
This construction does not require the input to be normalized (it will preform Trapezoidal Integration).
It does require the x array to be in ascending order.

# Examples
```julia-repl

julia> x,y = [1.0, 2.0, 3.0], [0.75, 0.5, 0.25]
([1.0, 2.0, 3.0], [0.75, 0.5, 0.25])

julia> LinearInterpolatedPDF(x,y)
LinearInterpolatedPDF{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Interpolations.Gridded{Interpolations.Linear},Tuple{Array{Float64,1}}},Interpolations.Gridded{Interpolations.Linear}}(
pdf_itp: 3-element extrapolate(interpolate((::Array{Float64,1},), ::Array{Float64,1}, Gridded(Interpolations.Linear())), Throw()) with element type Float64:
 0.75
 0.5
 0.25
cdf_itp: 3-element extrapolate(interpolate((::Array{Float64,1},), ::Array{Float64,1}, Gridded(Interpolations.Linear())), Throw()) with element type Float64:
 0.0
 0.625
 1.0
invcdf_itp: 3-element extrapolate(interpolate((::Array{Float64,1},), ::Array{Float64,1}, Gridded(Interpolations.Linear())), Throw()) with element type Float64:
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
struct LinearInterpolatedPDF{T,N,ITP,IT} <: ContinuousUnivariateDistribution
    pdf_itp::Extrapolation{T,N,ITP,IT,Flat{Nothing}}
    cdf_itp::Extrapolation{T,N,ITP,IT,Flat{Nothing}}
    invcdf_itp::Extrapolation{T,N,<:GriddedInterpolation,<:Gridded{Linear},Throw{Nothing}}

    function LinearInterpolatedPDF(x, y)
        area = integrate(x,y)
        p = y./area
        
        issorted(x) || error("x is not in ascending order")
        
        # pdf (zero outside of bounds)
        xmin = prevfloat(float(first(x)))
        xmax = nextfloat(float(last(x)))
        p_itp = LinearInterpolation([xmin; x; xmax], [0.0; p; 0.0],
                                    extrapolation_bc = Flat()
                                   )
        
        # cdf (0 below, 1 above)
        c = cumul_integrate(x,p)
        c_itp = LinearInterpolation([xmin; x; xmax], [0.0; c; 1.0],
                                    extrapolation_bc = Flat()
                                   )

        # invcdf (throw error out of bounds)
        i_itp = LinearInterpolation(c, x, extrapolation_bc=Throw())
        return LinearInterpolatedPDF(p_itp, c_itp, i_itp)
    end
end

pdf(d::LinearInterpolatedPDF, x::T) where {T<:Real} = d.pdf_itp(x)
cdf(d::LinearInterpolatedPDF, x::T) where {T<:Real} = d.cdf_itp(x)
quantile(d::LinearInterpolatedPDF, x::T) where {T<:Real} = d.invcdf_itp(x)

pdf(d::LinearInterpolatedPDF, x::AbstractVector{T}) where {T<:Real} = d.pdf_itp(x)
cdf(d::LinearInterpolatedPDF, x::AbstractVector{T}) where {T<:Real} = d.cdf_itp(x)
quantile(d::LinearInterpolatedPDF, x::AbstractVector{T}) where {T<:Real} = d.invcdf_itp(x)

using Random
import Base.rand
rand(rng::AbstractRNG, d::LinearInterpolatedPDF) = d.invcdf_itp(rand(rng))
