using Pkg
pkg"rm CardinalBSplines"
pkg"rm InfiniteVectors"
pkg"add https://github.com/vincentcp/InfiniteVectors.jl"
Pkg.develop(PackageSpec(path=splitdir(@__DIR__)[1]))
pkg"instantiate"
using Documenter, CardinalBSplines

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
        "Manual" => Any["man/evaluation.md",
            "man/integration.md",
            "man/filters.md"
            ]
        ],
    doctest=true
)

if "deploy" in ARGS && Sys.ARCH === :x86_64 && Sys.KERNEL === :Linux
    deploydocs(
        repo = "github.com/vincentcp/CardinalBSplines.jl.git",
    )
end
