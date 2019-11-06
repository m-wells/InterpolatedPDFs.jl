#=
    randpnt
    Copyright © 2019 m-wells <m-wells@pavilionhp>

    Distributed under terms of the MIT license.
=#

using Interpolations
using InterpolatedPDFs: randpnt, bin_mc, getknots

using Random: seed!

function fake_sinterp()
    seed!(42)
    sf1 = @addsupport (x,y) -> x+y<1.5
    sf2 = @addsupport (x,y) -> x+y>0.5
    knots = ([0; sort(rand(8)); 1],
             [0; sort(rand(7)); 1]
            )
    coeffs = rand(length.(knots)...)
    for cind in CartesianIndices(coeffs)
        if !(sf1(knots[cind]) && sf2(knots[cind]))
            coeffs[cind] = 0
        end
    end
    return SupportedInterp(knots, coeffs, [sf1, sf2])
end

function fake_interp()
    seed!(42)
    knots = ([0; sort(rand(8)); 1],
             [0; sort(rand(7)); 1]
            )
    coeffs = rand(length.(knots)...)
    return LinearInterpolation(knots, coeffs)
end

(interp::Interpolations.Extrapolation)(args::Tuple) = interp(args...)

function foo_mc(binds::CartesianIndices, sinterp::SupportedInterp)
    knots = getknots(sinterp)
    n = 1000

    num_current = 0
    sum_current = 0.0

    num_current += n
    pnt = randpnt(knots, binds)
    sum_current += sinterp(0.3,0.4)
end

foo(x,y) = x+y

function test()
    #x = [0; sort(rand(8)); 1]
    #inds = 2:3
    #randpnt(x, inds)
    #@time randpnt(x, inds)

    #binds = CartesianIndices((2:3,4:5))
    #sinterp = fake_sinterp()
    #seed!(42)
    #foo_mc(binds, sinterp)
    #@time foo_mc(binds, sinterp)
    
    #foo_mc(binds, sinterp)
    #@time foo_mc(binds, sinterp)

    #foo(1,2)
    #@time foo(1,2)
    #foo((1,2)...)
    #@time foo((1,2)...)

    #interp = fake_interp()
    #interp(0.4,0.3)
    #@time interp(0.4,0.3)
    #interp((0.4,0.3))
    #@time interp((0.4,0.3))
    
    #sinterp = fake_sinterp()
    #@code_warntype sinterp(0.4,0.3)
    #sinterp(0.4,0.3)
    #@time sinterp(0.4,0.3)
    #sinterp((0.4,0.3))
    #@time sinterp((0.4,0.3))


    sf1 = @addsupport (x,y) -> x+y<1.5
    @code_warntype sf1(0.3,0.2)
    @time sf1(0.3,0.2)

end
test()
