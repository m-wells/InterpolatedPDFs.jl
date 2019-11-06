using InterpolatedPDFs
using InterpolatedPDFs: randpnt, bin_mc, Δ, get_pmf, coefficients
using Test

import Random.seed!

function test()
    seed!(42)
    
    x = sort([0; rand(4); 1])
    y = sort([0; rand(3); 1])
    knots = (x,y)
    cinds = CartesianIndices((2:3,1:2))
    
    #@testset "randpnt" begin
    #    @test x[2] ≤ randpnt(x, 2:3) ≤ x[3]
    #    @test x[4] ≤ randpnt(x, 4:5) ≤ x[5]
    #    @test y[1] ≤ randpnt(y, 1:2) ≤ y[2]
    #
    #    pnt = randpnt(knots, cinds)
    #    @test x[2] ≤ pnt[1] ≤ x[3]
    #    @test y[1] ≤ pnt[2] ≤ y[2]
    #
    #    pnts = randpnts(knots, cinds, 5)
    #    for pnt in pnts
    #        @test x[2] ≤ pnt[1] ≤ x[3]
    #        @test y[1] ≤ pnt[2] ≤ y[2]
    #    end
    #end
    
    @testset "bin_mc" begin
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
        
        sinterp = SupportedInterp(knots, coeffs, [sf1, sf2])
        #@code_warntype sinterp(0.2,0.4)
        #sinterp(0.2,0.4)
        #@time sinterp(0.2,0.4)
        #@code_warntype sinterp(randpnt(knots, cinds))
        #sinterp(randpnt(knots, cinds))
        #@time sinterp(randpnt(knots, cinds))
        #@test isapprox(bin_mc(CartesianIndices((1:2,3:4)), sinterp), 0.103983747)
        #seed!(42)
        #@code_warntype bin_mc(CartesianIndices((1:2,3:4)), sinterp)
        #seed!(42)
        #bin_mc(CartesianIndices((1:2,3:4)), sinterp)
        seed!(42)
        val = bin_mc(CartesianIndices((1:2,3:4)), sinterp)
        @test isapprox(val, 0.10475865091044978, rtol=1e-3)
        #seed!(42)
        #@time bin_mc(CartesianIndices((1:2,3:4)), sinterp)
    
    end    

    @testset "Δ(knots)" begin
        knots = ([0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
                 [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
                )
        @test Δ(knots[1], 2:3) == (knots[1][3] - knots[1][2])
        @test Δ(knots[1], 2:4) == (knots[1][4] - knots[1][2])
        @test Δ(knots[2], 1:2) == (knots[2][2] - knots[2][1])
        @test Δ(knots[2], 6:7) == (knots[2][7] - knots[2][6])
        @test Δ(knots, CartesianIndices((2:3,1:2))) ==
              (knots[1][3] - knots[1][2])*(knots[2][2] - knots[2][1])
    
        knots = (range(0, stop=π, length=5),
                 range(1, stop=π, length=7)
                )
    
        @test Δ(knots[1], 2:3) == (knots[1][3] - knots[1][2])
        @test Δ(knots[1], 2:4) == (knots[1][4] - knots[1][2])
        @test Δ(knots[2], 1:2) == (knots[2][2] - knots[2][1])
        @test Δ(knots[2], 6:7) == (knots[2][7] - knots[2][6])
    end
    
    @testset "get_pmf" begin
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
        
        sinterp = SupportedInterp(knots, coeffs, [sf1, sf2])
    
        #@code_warntype get_pmf(sinterp)
        #@show get_pmf(sinterp)
        #@time get_pmf(sinterp)
        pmf = get_pmf(sinterp)
        total = sum(pmf)
    
        coeffs = coefficients(sinterp)
        coeffs ./= total
    
        pmf = get_pmf(sinterp)
        display(pmf)
        @test isapprox(sum(pmf), 1.0, rtol=1e-3)
    end
end
test()
