
# CardinalBSplines.jl Documentation

*A Julia package for B-splines operations*

For installation instructions, see [Installation](@ref).

For a  full description of the functionality use the manual:
```@contents
Pages = ["man/evaluation.md","man/integration.md","man/filters.md"]
```

## Installation

CardinalBSplines.jl is not added to the Julia package manager and depends on a unregistered package Sequences.jl.
The packages can easily be installed by cloning their git repository. From the Julia REPL, type `]` to enter Pkg mode and run

```julia
pkg> add https://github.com/vincentcp/Sequences.jl
pkg> add https://github.com/vincentcp/CardinalBSplines.jl
```

or in a file you could use

```julia
using Pkg
pkg"add https://github.com/vincentcp/Sequences.jl"
pkg"add https://github.com/vincentcp/CardinalBSplines.jl"
```
