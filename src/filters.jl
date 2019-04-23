"""
Discrete B spline signal.

Implements the signals of `B-Spline Signal Processing: Part I`\\
Unser et al., IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 41, NO. 2
"""
bsplinesignal(n::Int, ::Type{T}=Float64) where {T} =
    CompactInfiniteVector(evaluate_centered_BSpline.(Val{n}(), LinRange(convert(T, -n), convert(T, n), 2n+1), T), -n)
bsplinesignal(n::Int, m::Int, ::Type{T}=Float64) where {T} =
    CompactInfiniteVector(evaluate_centered_BSpline.(Val{n}(), LinRange(convert(T, -n), convert(T, n), 2n*m+1), T), -n*m)
periodicbsplinesignal(n::Int, m::Int, N::Int, ::Type{T}=Float64) where {T} =
    PeriodicInfiniteVector(bsplinesignal(n, m, T), m*N)

leastsquares_dualperiodicbsplinesignal(n::Int, m::Int, N::Int, ::Type{T}=Float64) where {T} =
    leastsquares_inv(periodicbsplinesignal(n, m, N, T), m)

# leastsquares_dualbsplinesignal(n::Int, m::Int, ::Type{T}=Float64; opts...) where {T} =
#     leastsquares_inv(bsplinesignal(n, m, T), m)

function compact_dualperiodicbsplinesignal(n::Int, m::Int, N::Int, ::Type{T}=Float64) where {T}
    b = bsplinesignal(n, m, T)
    PeriodicInfiniteVector(inv(b, m), N*m)
end

function periodicleastsquarescoefficients(n::Int, m::Int, N::Int, ::Type{T}=Float64) where {T}
    b = periodicbsplinesignal(n, m, N, T)
    inv(downsample(b*b, m))
end

# function leastsquarescoefficients(n::Int, m::Int, ::Type{T}=Float64; opts...) where {T}
#     p = periodicleastsquarescoefficients(n, m, 10n*m+1, T)
#     @warn "Not done"
#     CompactInfiniteVector(p[-5nm:5nm], -5nm)
# end


bn(n::Int, ::Type{T}=Float64) where T = bsplinesignal(n, T)
bnm(n::Int, m::Int, ::Type{T}=Float64) where T = bsplinesignal(n, m, T)
bnm(n::Int, m::Int, N::Int, ::Type{T}=Float64) where T = periodicbsplinesignal(n, m, N, T)
# b̃nm(n::Int, m::Int, ::Type{T}=Float64) where T = leastsquares_dualbsplinesignal(n, m, T)
b̃nm(n::Int, m::Int, N::Int, ::Type{T}=Float64) where T = leastsquares_dualperiodicbsplinesignal(n, m, N, T)
# snm(n::Int, m::Int, ::Type{T}=Float64; opts...) where T = leastsquarescoefficients(n, m, T; opts...)
snm(n::Int, m::Int, N::Int, ::Type{T}=Float64) where T = periodicleastsquarescoefficients(n, m, N, T)


IMPLEMENTED_FILTERS = []
for n in 0:5
    bn = Symbol(string("bn",n))
    @eval $bn = bn($n)
    @eval push!(IMPLEMENTED_FILTERS,$bn)
    @eval export $bn
    for m in 1:4
        bnm = Symbol(string("bnm",n, m))
        @eval $bnm = bnm($n,$m)
        @eval push!(IMPLEMENTED_FILTERS,$bnm)
        @eval export $bnm

        # bnm = Symbol(string("b̃nm",n, m))
        # @eval $bnm = b̃nm($n,$m)
        # @eval push!(IMPLEMENTED_FILTERS,$bnm)
        # @eval export $bnm
        #
        # bnm = Symbol(string("snm",n, m))
        # @eval $bnm = b̃nm($n,$m)
        # @eval push!(IMPLEMENTED_FILTERS,$bnm)
        # @eval export $bnm
    end
end
