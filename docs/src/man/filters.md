# Filters
The packages supports discrete B-splines as defined by

*B-Spline Signal Processing: Part I`Unser et al.,
IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 41, NO. 2*

$$b^n_m(k) = \beta^n(k/m)$$


```@setup filters
using PGFPlotsX, CardinalBSplines, LaTeXStrings
savefigs = (figname, obj) -> begin
    pgfsave(figname * ".pdf", obj)
    pgfsave(figname * ".tex", obj);
    pgfsave(figname * ".tikz", obj;include_preamble=false);
    run(`pdf2svg $(figname * ".pdf") $(figname * ".svg")`)
    return nothing
end
```

In the figure below we show two examples of such discrete B-splines, namely, $$b^1_3$$ and $$b^3_2$$

```@example filters
P = @pgf GroupPlot({group_style={group_size="2 by 1",},},  # hide
    {legend_cell_align="left",mark_options={fill_opacity=0.2}}, # hide
    PlotInc(bsplinesignal(1,3,Float64)), # hide
    LegendEntry(L"b^1_3"), # hide
    {legend_cell_align="left",mark_options={fill_opacity=0.2}}, # hide
    PlotInc(bsplinesignal(3,2,Float64)), # hide
    LegendEntry(L"b^3_2")) # hide
savefigs("discretespline", P) # hide
```

[\[.pdf\]](discretespline.pdf), [\[generated .tex\]](discretespline.tex), [\[generated .tikz\]](discretespline.tikz)

![](discretespline.svg)


For oversampled discrete B-splines ($$m>1$$) there exists compact duals, i.e., signals, $$g$$, such that

$$\langle \delta_{lm} * b^n_m, \delta_{km} * g\rangle=\delta_{l-k},\quad \forall l,k\in \mathbb Z.$$

Their support is restricted to $$k=-R,\dots,R$$. Below you see some compact dual signals for the discrete B-splines above. 

```@example filters
P = @pgf GroupPlot({group_style={group_size="2 by 1",},}, # hide
    {legend_cell_align="left",mark_options={fill_opacity=0.2}}, # hide
    PlotInc(inv(bsplinesignal(1,3,Float64),3)), # hide
    LegendEntry(L"R=3"), # hide
    PlotInc(inv(bsplinesignal(1,3,Float64),3,R=11)), # hide
    LegendEntry(L"R=5"), # hide
    PlotInc(inv(bsplinesignal(1,3,Float64),3,R=21)), # hide
    LegendEntry(L"R=10"), # hide
    {legend_cell_align="left",mark_options={fill_opacity=0.2}}, # hide
    PlotInc(inv(bsplinesignal(3,2,Float64),2)), # hide
    LegendEntry(L"R=4"), # hide
    PlotInc(inv(bsplinesignal(3,2,Float64),2,R=15)), # hide
    LegendEntry(L"R=7"), # hide
    PlotInc(inv(bsplinesignal(3,2,Float64),2,R=21)), # hide
    LegendEntry(L"R=10"),) # hide
savefigs("compact_dual", P) # hide
```

[\[.pdf\]](compact_dual.pdf), [\[generated .tex\]](compact_dual.tex), [\[generated .tikz\]](compact_dual.tikz)

![](compact_dual.svg)

```@autodocs
Modules = [CardinalBSplines]
Pages   = ["filters.jl"]
```
