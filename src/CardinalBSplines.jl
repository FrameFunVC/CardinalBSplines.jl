# CardinalBSplines.jl
module CardinalBSplines
using SpecialFunctions: gamma


include("splinetypes.jl")
include("splineevaluation.jl")
include("integration.jl")

using InfiniteVectors
export bn, bnm, snm, b̃nm, compact_dualperiodicbsplinesignal
include("filters.jl")

end # module
