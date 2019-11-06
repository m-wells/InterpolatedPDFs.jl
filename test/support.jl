using InterpolatedPDFs
using Base.Iterators
using Test
using Suppressor

import Random.seed!
seed!(42)

sf1 = @supportfunction (x,y) -> x+y<1.5 Tuple{Float64,Float64}
sf2 = @supportfunction (x,y) -> x+y>0.5 NTuple{2,Float64}
sf3 = @suppress @supportfunction (x,y) -> x+y>0.5
@testset "SupportFunction type" begin
    @test isa(sf1, InterpolatedPDFs.SupportFunction)
    @test isa(sf2, InterpolatedPDFs.SupportFunction)
    @test isa(sf3, InterpolatedPDFs.SupportFunction)
    @test_logs (:warn, r".*abstract.*") @supportfunction (x,y) -> x+y>0.5
end

knots = ([0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
         [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        )

@testset "SupportFunction testing across knots" begin
    for knot in product(knots...)
        if (sum(knot) < 1.5) && (sum(knot) > 0.5)
            @test (sf1(knot) && sf2(knot))
        else
            @test !(sf1(knot) && sf2(knot))
        end
    end
end

coeffs = rand(length.(knots)...)
for cind in CartesianIndices(coeffs)
    if !(sf1(knots[cind]) && sf2(knots[cind]))
        coeffs[cind] = 0
    end
end

sinterp = SupportedInterp(knots, coeffs, [sf1, sf2])
@testset "SupportInterp calling" begin
    @test iszero(sinterp(0.2,0.2))
    @test !iszero(sinterp(0.5,0.5))
    @test iszero(sinterp(0.8,0.8))
end

@testset "SupportInterp getknots/coefficients" begin
    @test InterpolatedPDFs.getknots(sinterp) === knots
    @test InterpolatedPDFs.coefficients(sinterp) === coeffs
end

@testset "SupportInterp indexing" begin
    cinds = CartesianIndices((2:3,4:5))
    @test sinterp[cinds] == coeffs[cinds]
end
