@testset "conditional sample 2d" begin
    val = fake_val(8,8)
    inds = (2,:)

    #@code_warntype sample(val, inds)
    @inferred sample(val)
    for _ in 1:10
        @test sample(val, inds)[1] == inds[1]
        @test 1 ≤ sample(val, inds)[2] ≤ size(val,2)
    end
    #tfunc = x -> @time sample(x)
    #@show tfunc(val)
end
