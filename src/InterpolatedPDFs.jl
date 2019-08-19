module InterpolatedPDFs

using Distributions
using Interpolations
using Interpolations: Extrapolation
using NumericalIntegration

export LinearInterpolatedPDF, fit_cpl, pdf, cdf, quantile

include("linear_1d.jl")

end # module
