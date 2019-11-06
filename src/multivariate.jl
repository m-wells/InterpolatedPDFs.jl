function randpnt(x::AbstractVector, inds::UnitRange{Int})
    i1,i2 = extrema(inds)
    return x[i1] + rand()*(x[i2] - x[i1])
end

function randpnt(knots::NTuple{N,AbstractVector},
                 cinds::CartesianIndices{N,NTuple{N,UnitRange{Int}}}
                ) where N
    return ntuple(i -> randpnt(knots[i], cinds.indices[i]), Val(N))
end

#function randpnts(knots::NTuple{N,AbstractVector},
#                  cinds::CartesianIndices{N,NTuple{N,UnitRange{Int}}},
#                  n::Int
#                 ) where N
#    retval = Vector{NTuple{N,Float64}}(undef, n)
#    for i in 1:n
#        retval[i] = randpnt(knots, cinds)
#    end
#    return retval
#end

function Δ(knots::AbstractRange{T}, inds::UnitRange{Int}) where T
    return length(knots) > 1 ? (length(inds)-1)*step(knots) : one(T)
end

function Δ(knots::AbstractVector{T}, inds::UnitRange{Int}) where T
    i1,i2 = extrema(inds)
    return length(knots) > 1 ? knots[i2]-knots[i1] : one(T)
end

function Δ(knots::NTuple{N,AbstractVector},
           cinds::CartesianIndices{N,NTuple{N,UnitRange{Int}}}
          ) where N
    return mapreduce(Δ, *, knots, cinds.indices)
end

#mutable struct 

"""
    bin_mc(binds::CartesianIndices, sinterp::SupportedInterp)
    
Use a Monte-Carlo scheme on a given bin of `sinterp` with corners located at `binds` to determine an average bin value.
This is used for bins that are on boundaries.
Does not apply `Δxs`.
"""
function bin_mc(binds::CartesianIndices{N,NTuple{N,UnitRange{Int}}},
                sinterp::SupportedInterp
               ) where N
    knots = getknots(sinterp)
    n = 2^8

    num_current = n
    sum_current = 0.0
    for i = 1:n
        sum_current += sinterp(randpnt(knots,binds))
    end

    num_prev = 0
    sum_prev = 0.0

    while !isapprox(sum_current/num_current, sum_prev/num_prev)
        num_prev = num_current
        sum_prev = sum_current

        num_current += n
        for i = 1:n
            sum_current += sinterp(randpnt(knots,binds))
        end
    end

    return sum_current/num_current
end

"""
    bin_avg(binds::CartesianIndices, sinterp::SupportedInterp)
    
Determine the average value for a bin.
This is used for bins that are not on boundaries.
Does not apply `Δxs`.
"""
bin_avg(binds::CartesianIndices, s::SupportedInterp) = sum(s[binds])/length(binds)

function get_pmf(sinterp::SupportedInterp{T,N,ITP,IT,BC}) where {T,N,ITP,IT,BC}
    binds = eachbin(sinterp)
    pmf = Array{Float64}(undef, size(binds))
    
    knots = getknots(sinterp)
    for (i,bind) in enumerate(binds)
        Δx = Δ(knots, bind)
        if isboundary(sinterp, i)
            pmf[i] = Δx*bin_mc(bind, sinterp)
        else
            pmf[i] = Δx*bin_avg(bind, sinterp)
        end
    end
    return pmf
end

struct MultivariatePDF{T,N,ITP,IT,BC,S}
    pdf::SupportedInterp{T,N,ITP,IT,BC,S}
    pmf::Array{T,N}

    function MultivariatePDF(pdf::SupportedInterp{T,N,ITP,IT,BC,S}
                            ) where {T,N,ITP,IT,BC,S}
        pmf = get_pmf(pdf)
        return new{T,N,ITP,IT,BC}(pdf, pmf)
    end
end
