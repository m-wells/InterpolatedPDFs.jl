"""
    Δ(knots::AbstractRange{T}, inds::UnitRange{Int}) where T

Efficient computation of `knots[last(inds)] - knots[first(inds)]` when `knots` is a range.
Returns `one(T)` when `length(inds)==1`.
"""
function Δ(knots::AbstractRange{T}, inds::UnitRange{Int}) where T
    length(inds) > 1 || begin
        length(inds) == 1 && return one(T)
        error("""
              inds is ill-formed
              inds = $(inds)
              """
             )
    end
    return (length(inds)-1)*step(knots)
end

"""
    Δ(knots::AbstractVector{T}, inds::UnitRange{Int}) where T

Efficient computation of `knots[last(inds)] - knots[first(inds)]` when `knots` is a vector.
Returns `one(T)` when `length(inds)==1`.
"""
function Δ(knots::AbstractVector{T}, inds::UnitRange{Int}) where T
    length(inds) > 1 || begin
        length(inds) == 1 && return one(T)
        error("""
              inds is ill-formed
              inds = $(inds)
              """
             )
    end
    i1,i2 = extrema(inds)
    @inbounds return knots[i2]-knots[i1]
end

"""
    Δ(knots::NTuple{N,AbstractVector}, inds::CartesianIndices{N})

Efficient computation of the product of `Δ(knots[i], inds.indices[i])` for `i=1:N`.
"""
function Δ(knots::NTuple{N,AbstractVector}, cinds::CartUnitInds{N}) where N
    return mapreduce(Δ, *, knots, cinds.indices)
end

"""
    bin_mc(binds::CartesianIndices, sinterp::SupportedInterp)
    
Use a Monte-Carlo scheme on a given bin of `sinterp` with corners located at `binds` to determine an average bin value.
This is used for bins that are on boundaries.
Does not apply `Δxs`.
"""
function bin_mc(binds::CartUnitInds, sinterp::SupportedInterp)
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
bin_avg(binds::CartUnitInds, s::SupportedInterp) = sum(s[binds])/length(binds)

function MultivariatePDF(s::SupportedInterp{N,ITP,IT,BC,S}
                        ) where {N,ITP,IT,BC,S}
    binds = eachbin(s)
    pmf = Array{Float64}(undef, size(binds))

    knots = getknots(s)
    for (i,bind) in enumerate(binds)
        Δx = Δ(knots, bind)
        if isboundary(s, i)
            pmf[i] = Δx*bin_mc(bind, s)
        else
            pmf[i] = Δx*bin_avg(bind, s)
        end
    end

    norm_const = sum(pmf)
    pmf./=norm_const
    coefs_pdf = coefficients(s)
    coefs_pdf./=norm_const

    return MultivariatePDF{N,ITP,IT,BC,S}(s, PMF(pmf))
end

function MultivariatePDF(knots::NTuple{N,AbstractVector}, coeffs::AbstractArray{T,N},
                         sfs = SupportFunction{Tuple{Nothing}}[]
                        ) where {N,T}
    sinterp = SupportedInterp(knots, coeffs, sfs)
    return MultivariatePDF(sinterp)
end


#function sample(m::MultivariatePDF)
#    knots = getknots(m)
#    bind = sample(m.pmf)
#    maxp = max(
#    while true
#        pnt = randpnt(knots, bind)
#        sinterp(pnt) < 
#    end
#end
