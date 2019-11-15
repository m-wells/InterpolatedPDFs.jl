@testset "sample $(length(n))d" for n in [(64,), (8,8), (4,4,4)]
    val = fake_val(n...)
    pmf = fake_pmf(n...)
    sinterp = fake_sinterp()

    #@show sample(val)
    #@code_warntype sample(val)
    @inferred sample(val)
    @test sample(val) in indices(val)
    #tfunc = x -> @time sample(x)
    #tfunc(val)

    #@show samples(val, 100)
    #@code_warntype samples(val, 100)
    @inferred samples(val, 100)
    @test all(x -> x in indices(val), samples(val,100))
    #tfunc = x -> @time samples(x,100)
    #tfunc(val)

    #@show sample(pmf)
    #@code_warntype sample(pmf)
    @inferred sample(pmf)
    @test sample(pmf) in indices(pmf)
    #tfunc = x -> @time sample(x)
    #tfunc(pmf)

    #@show samples(pmf, 100)
    #@code_warntype samples(pmf, 100)
    @inferred samples(pmf, 100)
    @test all(x -> x in indices(pmf), samples(pmf,100))
    #tfunc = x -> @time samples(x,100)
    #tfunc(pmf)
end

@testset "sample SupportedInterp" begin
    sinterp = fake_sinterp()
    cinds = CartesianIndices((2:3,1:2))
    @test_throws DomainError sample(sinterp, cinds)

    i = 4:5
    j = 3:4
    cinds = CartesianIndices((i,j))
    x,y = getknots(sinterp)
    pnt1 = (x[first(i)], y[first(j)])
    pnt2 = (x[last(i)], y[last(j)])
    for _ = 1:10
        @test is_between(sample(sinterp, cinds), pnt1, pnt2)
    end
end

@testset "sample MultivariatePDF" begin
    mpdf = fake_mpdf()
    cinds = CartesianIndices((2:3,1:2))
    @test_throws DomainError sample(sinterp, cinds)

    i = 4:5
    j = 3:4
    cinds = CartesianIndices((i,j))
    x,y = getknots(sinterp)
    pnt1 = (x[first(i)], y[first(j)])
    pnt2 = (x[last(i)], y[last(j)])
    for _ = 1:10
        @test is_between(sample(sinterp, cinds), pnt1, pnt2)
    end
end
