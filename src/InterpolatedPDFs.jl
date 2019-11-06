module InterpolatedPDFs

using FunctionWrappers: FunctionWrapper
#using Distributions
#using NumericalIntegration
using MacroTools
using Base.Iterators
using Interpolations
using Interpolations: Extrapolation
import Interpolations: coefficients, getknots

coefficients(etp::Extrapolation) = coefficients(etp.itp)

export @supportfunction, SupportedInterp
#export LinearInterpolatedPDF, pdf, cdf, quantile

#include("./linear_1d.jl")

include("./support.jl")
include("./cartesian_binning.jl")
include("./multivariate.jl")

end # module
