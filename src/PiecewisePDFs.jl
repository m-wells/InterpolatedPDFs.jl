module PiecewisePDFs

using Distributions
using Interpolations
using Interpolations: Extrapolation
using NumericalIntegration

export ContinuousPiecewiseLinear, fit_cpl, pdf, cdf, quantile

include("piecewise_1d.jl")

end # module
