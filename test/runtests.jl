using InterpolatedPDFs
using Random
using Test

Random.seed!(1234)

@testset "linear_1d" begin

    @testset "UnitRange knots" begin
        x = 0:10
        d1 = fit_pdf(x, 10.0.*rand(10))
        p = pdf.(d1, x)
        d2 = LinearInterpolatedPDF(x, p)

        for d in [d1, d2]
            @test getknots(d) == x

            @test iszero(cdf(d, first(x)))
            @test isone(cdf(d, last(x)))

            @test quantile(d, 0) == first(x)
            @test quantile(d, 1) == last(x)
        end
    end

    @testset "StepRangeLen knots" begin
        x = range(0, stop=5, length=10)
        d1 = fit_pdf(x, 5.0.*rand(10))
        p = pdf.(d1, x)
        d2 = LinearInterpolatedPDF(x, p)

        for d in [d1, d2]
            @test getknots(d) == x

            @test iszero(cdf(d, first(x)))
            @test isone(cdf(d, last(x)))

            @test quantile(d, 0) == first(x)
            @test quantile(d, 1) == last(x)
        end
    end

    @testset "Vector knots" begin
        x = [0.0, 0.5, 1.0, 1.5, 2.0]
        d1 = fit_pdf(x, 2.0.*rand(10))
        p = pdf.(d1, x)
        d2 = LinearInterpolatedPDF(x, p)

        for d in [d1, d2]
            @test getknots(d) == x

            @test iszero(cdf(d, first(x)))
            @test isone(cdf(d, last(x)))

            @test quantile(d, 0) == first(x)
            @test quantile(d, 1) == last(x)

            @test iszero(cdf(d, first(x)))
            @test isone(cdf(d, last(x)))

            @test quantile(d, 0) == first(x)
            @test quantile(d, 1) == last(x)
        end
    end

    @testset "pdf cdf quantile check" begin
        x = [0.0, 0.5, 1.0, 1.5, 2.0]
        d1 = fit_pdf(x, 2.0.*rand(10))
        p = pdf.(d1, x)
        d2 = LinearInterpolatedPDF(x, p)

        for d in [d1, d2]
            y = 0.2
            @test isa(pdf(d,y), Float64)

            @test isa(cdf(d,y), Float64)

            @test isa(quantile(d,y), Float64)

            n = 10
            y = pdf.(Ref(d), 2 .* rand(n))
            @test isa(y, Vector{Float64})
            @test length(y) == n

            y = cdf.(Ref(d), 2 .* rand(n))
            @test isa(y, Vector{Float64})
            @test length(y) == n

            y = quantile.(Ref(d), rand(n))
            @test isa(y, Vector{Float64})
            @test length(y) == n
        end
    end

    # TODO: better error checks
    @testset "error check" begin
        @test_throws ErrorException fit_pdf([0], rand(10))

        @test_throws ErrorException fit_pdf([1,2,3], 3.0.*rand(10))

        @test_throws ErrorException fit_pdf([0,1,2], 3.0.*rand(10))

        @test_throws ErrorException fit_pdf([0,2,1], 2.0.*rand(10))
    end

    @testset "accuracy check" begin
        x = range(0, stop=π/2, length=20)
        d1 = fit_pdf(x, acos.(rand(100000)))
        p = pdf.(d1, x)
        d2 = LinearInterpolatedPDF(x, p)
        xmid = InterpolatedPDFs.midpoints(x)

        for d in [d1, d2]
            @test maximum(abs.(pdf.(Ref(d), xmid) .- sin.(xmid))) < 0.02
        end
    end

    @testset "sampling" begin
        x = range(0, stop=π/2, length=10)
        d1 = fit_pdf(x, acos.(rand(100)))
        p = pdf.(d1, x)
        d2 = LinearInterpolatedPDF(x, p)

        for d in [d1, d2]
            s = rand(d,100)
            @test minimum(s) ≥ first(x)
            @test maximum(s) ≤ last(x)
            @test isa(s,Vector{Float64})

            s = rand(d,10,10)
            @test isa(s,Matrix{Float64})

            s = rand(d,5,4,3)
            @test isa(s,Array{Float64,3})
        end
    end
end
