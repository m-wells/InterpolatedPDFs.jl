using InterpolatedPDFs
using InterpolatedPDFs: randpnt, bin_mc, Δ, coefficients, PMF
import Random.seed!
using Test

indices(x::AbstractVector) = eachindex(x)
indices(x::AbstractArray) = CartesianIndices(x)
indices(pmf::PMF) = indices(pmf.vals)

function fake_knots(n::NTuple{N,Int}) where N
    seed!(42)
    return ntuple(i->sort([0; rand(n[i]-2); 1]), Val(N))
end
fake_knots(n...) = fake_knots(n)

function fake_val(n...)
    seed!(42)
    return rand(n...)
end

function fake_pmf(n...)
    val = fake_val(n...)
    val./=sum(val)
    return PMF(val)
end

function fake_sinterp()
    sf1 = @supportfunction (x,y) -> x+y<1.5 Float64
    sf2 = @supportfunction (x,y) -> x+y>0.5 Float64
    knots = ([0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
             [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
            )
    seed!(42)
    coeffs = rand(length.(knots)...)
    for cind in CartesianIndices(coeffs)
        if !(sf1(knots[cind]) && sf2(knots[cind]))
            coeffs[cind] = 0
        end
    end
    
    return SupportedInterp(knots, coeffs, [sf1, sf2])
end

function fake_mpdf()
    return MultivariatePDF(fake_sinterp())
end

is_between(x::T, x1::T, x2::T) where T<:Real = x1 < x2 ? x1 ≤ x ≤ x2 : x2 ≤ x ≤ x1

function is_between(x::T, x1::T, x2::T) where {N,T<:NTuple{N,Real}}
    for i = 1:N
        @inbounds is_between(x[i], x1[i], x2[i]) || return false
    end
    return true
end



#include("./cartesian_binning.jl")
#include("./support.jl")
#include("./multivariate.jl")
#include("./sample.jl")
include("./conditional.jl")
