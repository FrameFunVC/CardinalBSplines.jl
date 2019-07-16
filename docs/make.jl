using Pkg
Pkg.develop(PackageSpec(path=splitdir(@__DIR__)[1]))
using Documenter, CardinalBSplines, PGFPlotsX, LaTeXStrings, DocumentPGFPlots

imgdir = normpath(joinpath(@__DIR__(),"src","man","figs"))

P = @pgf GroupPlot({group_style={group_size="2 by 1",},},  # hide
                {legend_cell_align="left",mark_options={fill_opacity=0.2}}, # hide
                PlotInc(bsplinesignal(1,3,Float64)), # hide
                LegendEntry(L"b^1_3"), # hide
                {legend_cell_align="left",mark_options={fill_opacity=0.2}}, # hide
                PlotInc(bsplinesignal(3,2,Float64)), # hide
                LegendEntry(L"b^3_2")) # hide
DocumentPGFPlots.savefigs(joinpath(imgdir,"discretespline"), P)
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
DocumentPGFPlots.savefigs(joinpath(imgdir,"compact_dual"), P)

const render_pdf = "pdf" in ARGS
let r = r"buildroot=(.+)", i = findfirst(x -> occursin(r, x), ARGS)
    global const buildroot = i === nothing ? (@__DIR__) : first(match(r, ARGS[i]).captures)
end

const format = if render_pdf
    LaTeX(
        platform = "texplatform=docker" in ARGS ? "docker" : "native"
    )
else
    Documenter.HTML(
        prettyurls = ("deploy" in ARGS),
    )
end

makedocs(sitename="CardinalBSplines.jl",
    modules = [CardinalBSplines],
    authors = "vincentcp",
    format = format,
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "man/evaluation.md",
            "man/integration.md",
            "man/filters.md",
            ]
        ],
    doctest= ("deploy" in ARGS) ? true : :fix
)

if "deploy" in ARGS && Sys.ARCH === :x86_64 && Sys.KERNEL === :Linux
    deploydocs(
        repo = "github.com/FrameFunVC/CardinalBSplines.jl.git",
    )
end
