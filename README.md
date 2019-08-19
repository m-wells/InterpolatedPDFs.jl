# PiecewisePDFs.jl

Simple extension of [Distributions.jl](https://github.com/JuliaStats/Distributions.jl) providing support for piecewise pdfs.
Currently only one type is implemented

```
ContinuousPiecewiseLinear{T,1} <: ContinuousUnivariateDistribution
```

A continuous univariate linear piecewise distribution.
The pdf, cdf, and inverse cdf are interpolated using [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl).

# Examples
The easiest way to create a distribution is to use `fit_cpl`
```
julia> x = range(0, pi/2, length=10)
0.0:0.17453292519943295:1.5707963267948966

julia> s = acos.(rand(1000));

julia> d = fit_cpl(x,s)
ContinuousPiecewiseLinear{Float64,1}(
pdf_itp: 10-element extrapolate(scale(interpolate(::Array{Float64,1}, BSpline(Interpolations.Linear())), (0.0:0.17453292519943295:1.5707963267948966,)), Throw()) with element type Float64:
 0.0575436610456399
 0.18991497776431412
 0.44360944378499123
 0.5467161666950493
 0.7497086195023982
 0.6657386372811903
 0.8012734435737708
 0.9451934695901671
 1.0121028020550467
 0.6930971210769691
cdf_itp: 10-element extrapolate(scale(interpolate(::Array{Float64,1}, BSpline(Interpolations.Linear())), (0.0:0.17453292519943295:1.5707963267948966,)), Throw()) with element type Float64:
 0.0
 0.0215948400486856
 0.07688027528782507
 0.1633024881363229
 0.27643689325436793
 0.39995796835034425
 0.5279789232376059
 0.6803869127968689
 0.8511932346829605
 1.0
invcdf_itp: 10-element extrapolate(interpolate((::Array{Float64,1},), ::Array{Float64,1}, Gridded(Interpolations.Linear())), Throw()) with element type Float64:
 0.0
 0.17453292519943295
 0.3490658503988659
 0.5235987755982988
 0.6981317007977318
 0.8726646259971648
 1.0471975511965976
 1.2217304763960306
 1.3962634015954636
 1.5707963267948966
)
```

After fitting the distribution you can do useful things like
```
julia> pdf(d,1)
0.7646218435870893

julia> cdf(d,0.5)
0.1516172522108395

julia> quantile(d,0.9)
1.4535080268795113

julia> rand(d,10)
10-element Array{Float64,1}:
 0.27565417806686643
 1.074337923663701
 1.237530643864552
 0.4744230962935516
 1.18776692814955
 0.8436400094154567
 1.0835325983972564
 1.1413257453616537
 0.8701141622223004
 1.1702951450424084
```
