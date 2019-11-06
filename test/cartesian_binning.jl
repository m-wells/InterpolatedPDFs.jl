using InterpolatedPDFs: EachBin, eachbin 
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
