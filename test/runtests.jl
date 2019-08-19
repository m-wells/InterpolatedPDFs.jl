using InterpolatedPDFs
using Test

import Random.seed!
seed!(1234)

@testset "linear_1d" begin
    x = range(0.0, π/2, length=10)
    d = fit_cpl(x, acos.(rand(10000)))
    @test isa(d, LinearInterpolatedPDF{Float64,1})

    r = maximum(abs.(pdf(d,x) .- sin.(x)))
    @test r < 0.05

    @test iszero(cdf(d, first(x)))
    @test isone(cdf(d, last(x)))

    @test quantile(d, 0) == first(x)
    @test quantile(d, 1) == last(x)

    s = rand(d)
    @test s ≥ first(x)
    @test s ≤ last(x)
    @test isa(s,Float64)

    s = rand(d,100)
    @test minimum(s) ≥ first(x)
    @test maximum(s) ≤ last(x)
    @test isa(s,Vector{Float64})

    s = rand(d,10,10)
    @test isa(s,Matrix{Float64})

    s = rand(d,5,4,3)
    @test isa(s,Array{Float64,3})
end
