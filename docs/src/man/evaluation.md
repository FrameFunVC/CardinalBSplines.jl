```@setup manual
using Plots, CardinalBSplines
```
# B-spline evaluation
Spline evalution can be done in two ways. The first one uses the spline constructors in [B-spline constructors](@ref).

```@example manual
f = BSpline(3)
x = LinRange(-1,5,10)
f.(x)
```

The second one is uses the internal function in [Internal evaluation functions](@ref).

```@example manual
CardinalBSplines.evaluate_Bspline.(Val{3}(), x, Float64)
```

```@meta
DocTestFilters = r"[0-9\.]+ seconds \(.*\)"
```

No allocation happens during evaluation
```@jldoctest manual
x = LinRange(-1,5,10000);
s = Array{Float64}(undef, 10000);
s .= f.(x);
@time s .= f.(x);
  0.000165 seconds (5 allocations: 208 bytes)
```
```@meta
DocTestFilters = nothing
```
## Examples
### B-splines
```@example manual
x = LinRange(-1,11,500) # hide
plot(legend=false)
[plot!(x,BSpline(degree).(x)) for degree in 0:10]
Plots.savefig("bspline-plot.svg"); nothing # hide
```
![](bspline-plot.svg)
### Centered B-splines
```@example manual
x = LinRange(-6,6,500) # hide
plot(legend=false)
[plot!(x,CenteredBSpline(degree).(x)) for degree in 0:10]
Plots.savefig("cbspline-plot.svg"); nothing # hide
```
![](cbspline-plot.svg)
### Periodic B-splines
```@example manual
x = LinRange(-1,11,500) # hide
plot(legend=false)
[plot!(x,PeriodicBSpline(degree, 3).(x)) for degree in 0:10]
Plots.savefig("pbspline-plot.svg"); nothing # hide
```
![](pbspline-plot.svg)


### Derivatives of B-Splines
```@example manual
x = LinRange(-1,11,500) # hide
plot(legend=false)
[plot!(x,BSplineDiff(4,diff).(x)) for diff in 0:4]
Plots.savefig("dbspline-plot.svg"); nothing # hide
```
![](dbspline-plot.svg)




## B-spline constructors
```@autodocs
Modules = [CardinalBSplines]
Pages   = ["splinetypes.jl"]
```

## Internal evaluation functions
```@autodocs
Modules = [CardinalBSplines]
Pages   = ["splineevaluation.jl"]
```
