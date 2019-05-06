
# CardinalBSplines.jl Documentation

*A Julia package for B-splines operations*

For installation instructions, see [Installation](@ref).

For a  full description of the functionality use the manual:
```@contents
Pages = ["man/evaluation.md","man/integration.md","man/filters.md"]
```

## Installation

CardinalBSplines.jl is not added to the Julia General registry and depends on a unregistered package InfiniteVectors.jl.

### Recomanded
For Julia 1.1 or higher, you can add the FrameFun registry and than add CardinalBSplines.
From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> add https://github.com/vincentcp/FrameFunRegistry
pkg> add CardinalBSplines
```

### Legacy
In Julia 1.0, the packages can be installed by cloning their git repository. From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> add https://github.com/vincentcp/InfiniteVectors.jl
pkg> add https://github.com/vincentcp/CardinalBSplines.jl
```

or in a file you could use

```julia
using Pkg
pkg"add https://github.com/vincentcp/InfiniteVectors.jl"
pkg"add https://github.com/vincentcp/CardinalBSplines.jl"
```
