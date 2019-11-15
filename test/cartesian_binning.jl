using InterpolatedPDFs: EachBin, eachbin, randpnt, getknots
seed!(42)
knots = ([0; rand(3); 1],
         [0; rand(2); 1],
         [0; rand(1); 1]
        )
x = rand(5,4,3)
bins = eachbin(x)

@testset "EachBin iteration object" begin
    @test isa(bins, EachBin)
end

@testset "knot indexing" begin
    for inds in bins
        @test isa(inds, CartesianIndices)
        @test isa(knots[inds], Array{NTuple{3,Float64},3})
    end
end

@testset "randpnt" begin
    for _ = 1:10
        x1,x2 = 0.2, 1.7
        @test is_between(randpnt(x1,x2), x1, x2)
        x1,x2 = -1.7, -0.2
        @test is_between(randpnt(x1,x2), x1, x2)
        x1,x2 = -1.7, 0.2
        @test is_between(randpnt(x1,x2), x1, x2)

        x1,x2 = 1.7, 0.2
        @test is_between(randpnt(x1,x2), x1, x2)
        x1,x2 = -0.2, -1.7
        @test is_between(randpnt(x1,x2), x1, x2)
        x1,x2 = 1.7, -0.2
        @test is_between(randpnt(x1,x2), x1, x2)
    end

    x = rand(10).*5
    sort!(x)
    i = 2:3
    x1,x2 = x[i]
    for _ = 1:10
        @test is_between(randpnt(x,i), x1, x2)
    end
    for _ = 1:10
        @test !is_between(randpnt(x,i), x[4], x[5])
    end

    y = rand(9).*10
    sort!(y)
    j = 3:4
    c = CartesianIndices((i,j))
    pnt1 = (x[first(i)], y[first(j)])
    pnt2 = (x[last(i)], y[last(j)])
    for _ = 1:10
        @test is_between(randpnt((x,y), c), pnt1, pnt2)
    end

    sinterp = fake_sinterp()
    x,y = getknots(sinterp)
    pnt1 = (x[first(i)], y[first(j)])
    pnt2 = (x[last(i)], y[last(j)])
    for _ = 1:10
        @test is_between(randpnt(sinterp, c), pnt1, pnt2)
    end
end
