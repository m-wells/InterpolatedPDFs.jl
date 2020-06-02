module InterpolatedPDFs

using Distributions
using NumericalIntegration
using StatsBase: midpoints

using Interpolations
using Interpolations: getknots

using Random: AbstractRNG

export LinearInterpolatedPDF, fit_cpl, pdf, cdf, quantile, getknots

include("linear_1d.jl")

end # module
