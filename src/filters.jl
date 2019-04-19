"""
Discrete B spline signal.

Implements the signals of `B-Spline Signal Processing: Part I`\\
Unser et al., IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 41, NO. 2
"""
bsplinesignal(n::Int, ::Type{T}=Float64) where {T} =
    CompactInfiniteVector(evaluate_centered_BSpline.(Val{n}(), LinRange(convert(T, -n), convert(T, n), 2n+1), T), -n)
bsplinesignal(n::Int, m::Int, ::Type{T}=Float64) where {T} =
    CompactInfiniteVector(evaluate_centered_BSpline.(Val{n}(), LinRange(convert(T, -n), convert(T, n), 2n*m+1), T), -n*m)


IMPLEMENTED_FILTERS = []
for n in 0:5
    bn = Symbol(string("bn",n))
    # @eval $bn() = bsplinesignal($n)
    @eval $bn = bsplinesignal($n)
    @eval push!(IMPLEMENTED_FILTERS,$bn)
    @eval export $bn
    for m in 1:4
        bnm = Symbol(string("bnm",n, m))
        # @eval $bnm() = bsplinesignal($n,$m)
        @eval $bnm = bsplinesignal($n,$m)
        @eval push!(IMPLEMENTED_FILTERS,$bnm)
        @eval export $bnm
    end
end
