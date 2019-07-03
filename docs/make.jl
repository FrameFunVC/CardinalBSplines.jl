using Pkg
Pkg.develop(PackageSpec(path=splitdir(@__DIR__)[1]))
using Documenter, CardinalBSplines, PGFPlotsX, LaTeXStrings

imgdir = normpath(joinpath(@__DIR__(),"src","man","figs"))


savefigs = (figname, obj) -> begin
    pgfsave(figname * ".tex", obj);
    pgfsave(figname * ".tikz", obj;include_preamble=false);
    if "docker" in ARGS
        compile_tex(figname * ".tex", "docker")
        pdf2svg(figname, "docker")
    elseif "native" in ARGS
        compile_tex(figname * ".tex", "native")
        pdf2svg(figname, "native")
    else
        compile_tex(figname * ".tex")
        pdf2svg(figname)
    end
    return nothing
end

figname=joinpath(imgdir,"discretespline")

function piperun(cmd)
    verbose = "verbose" in ARGS || get(ENV, "DOCUMENTER_VERBOSE", "false") == "true"
    run(pipeline(cmd, stdout = verbose ? stdout : "LaTeXWriter.stdout",
                      stderr = verbose ? stderr : "LaTeXWriter.stderr"))
end

function compile_tex(texfile::String, platform = (Sys.which("latexmk") === nothing) ? "docker" : "native")
    thisdir = pwd()
    if platform == "native"
        Sys.which("latexmk") === nothing && (@error "LaTeXWriter: latexmk command not found."; return false)
        @info "LaTeXWriter: using latexmk to compile tex."
        try
            p = splitdir(figname)[1]=="" ? pwd() : splitdir(figname)[1]
            cd(p)
            piperun(`latexmk -f -interaction=nonstopmode -view=none -lualatex -shell-escape $(Base.basename(texfile))`)
            return true
        catch err
            logs = cp(pwd(), mktempdir(); force=true)
            @error "LaTeXWriter: failed to compile tex with latexmk. "
            return false
        end
    else
        Sys.which("docker") === nothing && (@error "LaTeXWriter: docker command not found."; return false)
        @info "LaTeXWriter: using docker to compile tex."
        script = """
            mkdir /home/zeptodoctor/build
            cd /home/zeptodoctor/build
            cp -r /mnt/. .
            latexmk -f -interaction=nonstopmode -view=none -lualatex -shell-escape $(Base.basename(texfile))
        """
        try
            p = splitdir(figname)[1]=="" ? pwd() : splitdir(figname)[1]
            @show p
            piperun(`docker run -itd -u zeptodoctor --name latex-container -v $(p):/mnt/ --rm juliadocs/documenter-latex:0.1`)
            piperun(`docker exec -u zeptodoctor latex-container bash -c $(script)`)
            piperun(`docker cp latex-container:/home/zeptodoctor/build/. $(p)`)
            return true
        catch err
            logs = cp(pwd(), mktempdir(); force=true)
            @error "LaTeXWriter: failed to compile tex with docker. Look for LaTeXWriter.stderr, LaTeXWriter.stdout"
            return false
        finally
            cd(thisdir)
            try; piperun(`docker stop latex-container`); catch; end
        end
    end
end

function pdf2svg(figname::String, platform = (Sys.which("pdf2svg") === nothing ) ? "docker" : "native")
    if platform=="native"
        Sys.which("pdf2svg") === nothing && (@error "pdf2svg."; return false)
        @info "using pdf2svg to transform pdf to svg."
        try
            run(`pdf2svg $(figname * ".pdf") $(figname * ".svg")`)
            return true
        catch err
            logs = cp(pwd(), mktempdir(); force=true)
            @error "failed to create svg with pdf2svg. "
            return false
        end
    else
        Sys.which("docker") === nothing && (@error "docker command not found."; return false)
        @info "using docker to create svg."
        try
            p = splitdir(figname)[1]=="" ? pwd() : splitdir(figname)[1]
            piperun(`docker run -itd --name svg-container -v $(p):/mnt/  --rm vincentcoppe/pdf2svg pdf2svg /mnt/$(Base.basename(figname)).pdf /mnt/$(Base.basename(figname)).svg`)
            return true
        catch err
            @error "failed to create svg with docker."
            return false
        finally
        end
    end
end




P = @pgf GroupPlot({group_style={group_size="2 by 1",},},  # hide
                {legend_cell_align="left",mark_options={fill_opacity=0.2}}, # hide
                PlotInc(bsplinesignal(1,3,Float64)), # hide
                LegendEntry(L"b^1_3"), # hide
                {legend_cell_align="left",mark_options={fill_opacity=0.2}}, # hide
                PlotInc(bsplinesignal(3,2,Float64)), # hide
                LegendEntry(L"b^3_2")) # hide
savefigs(joinpath(imgdir,"discretespline"), P)
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
savefigs(joinpath(imgdir,"compact_dual"), P)





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
        repo = "github.com/vincentcp/CardinalBSplines.jl.git",
    )
end
