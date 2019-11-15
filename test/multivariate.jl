n = (10,9)
@testset "randpnt" begin
    knots = fake_knots(n)
    x,y = knots

    @test x[2] ≤ randpnt(x, 2:3) ≤ x[3]
    @test x[4] ≤ randpnt(x, 4:5) ≤ x[5]
    @test y[1] ≤ randpnt(y, 1:2) ≤ y[2]

    cinds = CartesianIndices((2:3,1:2))
    pnt = randpnt(knots, cinds)
    @test x[2] ≤ pnt[1] ≤ x[3]
    @test y[1] ≤ pnt[2] ≤ y[2]
end

@testset "bin_mc" begin
    sinterp = fake_sinterp()
    seed!(42)
    val = bin_mc(CartesianIndices((1:2,3:4)), sinterp)
    @test isapprox(val, 0.10475579070838298, rtol=1e-3)
                        
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
    sinterp = fake_sinterp()
    pmf = get_pmf(sinterp)
    @test isapprox(sum(pmf), 1)
    pmf2 = get_pmf(sinterp)
    for i in eachindex(pmf)
        @test isapprox(pmf[i], pmf2[i], rtol=1e-2)
    end
end

@testset "MultivariatePDF" begin
    mpdf = fake_mpdf()
    @test isa(mpdf, MultivariatePDF)
end
