export bsplinesignal, periodicbsplinesignal
"""
    bsplinesignal(n::Int, [, m::Int] ::Type{T}=Float64)

Discrete B-spline signal of degree `n` and oversampling `m`.

Implements the signals of `B-Spline Signal Processing: Part I`\\
Unser et al., IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 41, NO. 2
"""
bsplinesignal(n::Int, ::Type{T}=Float64) where {T} =
    CompactInfiniteVector(evaluate_centered_BSpline.((Val{n}()),
        (LinRange(convert(T, -((n+2)>>1)), convert(T, ((n+2)>>1)-1), ((n+2)>>1)<<1)), T), -((n+2)>>1))
function bsplinesignal(n::Int, m::Int, ::Type{T}=Float64) where {T}
    K = (m*(n+1)+1)>>1
    if n==0
        CompactInfiniteVector(evaluate_centered_BSpline.(Val{n}(), LinRange(convert(T, -K)/m, convert(T, K-1)/m, K<<1), T), -K)
    else
        CompactInfiniteVector(evaluate_centered_BSpline.(Val{n}(), LinRange(convert(T, -K+1)/m, convert(T, K-1)/m, (K<<1)-1), T), -K+1)
    end
end

"""
    periodicbsplinesignal(n::Int, m::Int, N::Int, ::Type{T}=Float64)

Discrete B-spline signal of degree `n`, oversampling `m` and period `mN`.

See also `bsplinesignal`
"""
periodicbsplinesignal(n::Int, m::Int, N::Int, ::Type{T}=Float64) where {T} =
    PeriodicInfiniteVector(bsplinesignal(n, m, T), m*N)

leastsquares_dualperiodicbsplinesignal(n::Int, m::Int, N::Int, ::Type{T}=Float64) where {T} =
    leastsquares_inv(periodicbsplinesignal(n, m, N, T), m)

# leastsquares_dualbsplinesignal(n::Int, m::Int, ::Type{T}=Float64; opts...) where {T} =
#     leastsquares_inv(bsplinesignal(n, m, T), m)

function minimalK(p,q)
    max(0,iseven(p) ?
    ceil(Int, (p+1)/2*q/(q-1)-(q+1)/(q-1)) +1 :
    round((p+1)/2*q/(q-1)-(q+1)/(q-1)) ≈ (p+1)/2*q/(q-1)-(q+1)/(q-1) ?
        round(Int,(p+1)/2*q/(q-1)-(q+1)/(q-1)) + 1 :
        ceil(Int, (p+1)/2*q/(q-1)-(q+1)/(q-1)))
end

function compact_dualperiodicbsplinesignal(n::Int, m::Int, N::Int, ::Type{T}=Float64) where {T}
    b = bsplinesignal(n, m, T)
    PeriodicInfiniteVector(inv(b, m;K=minimalK(n,m)), N*m)
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

"""
    bn(n::Int, ::Type{T}=Float64)

Discrete B-spline signal of degree `n` and period `mN`.

Implements the signals of `B-Spline Signal Processing: Part I`\\
Unser et al., IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 41, NO. 2
"""
bn(n::Int, ::Type{T}=Float64) where T = bsplinesignal(n, T)
"""
    bnm(n::Int, m::Int, [N::Int, ]::Type{T}=Float64)

Discrete B-spline signal of degree `n` and oversampling `m`.
If 'N' is specified the signal is periodic with period `Nm`.

Implements the signals of `B-Spline Signal Processing: Part I`\\
Unser et al., IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 41, NO. 2
"""
bnm(n::Int, m::Int, ::Type{T}=Float64) where T = bsplinesignal(n, m, T)
bnm(n::Int, m::Int, N::Int, ::Type{T}=Float64) where T = periodicbsplinesignal(n, m, N, T)
# b̃nm(n::Int, m::Int, ::Type{T}=Float64) where T = leastsquares_dualbsplinesignal(n, m, T)
"""
    b̃nm(n::Int, m::Int, N::Int, ::Type{T}=Float64)

Least squares dual of `bnm`:

    \$b̃nm(k)= \\left[\\left([bnm*bnm]_{↓m}\\right)^{-1}\\right]_{↑m}*bnm(k)\$.
"""
b̃nm(n::Int, m::Int, N::Int, ::Type{T}=Float64) where T = leastsquares_dualperiodicbsplinesignal(n, m, N, T)
# snm(n::Int, m::Int, ::Type{T}=Float64; opts...) where T = leastsquarescoefficients(n, m, T; opts...)
"""
    snm(n::Int, m::Int, N::Int, ::Type{T}=Float64)

The coefficients of `b̃nm` in `bnm`, i.e.,

    \$\\tilde bnm(k) = [snm]_{\\uparrow m}*bnm(k)\$
"""
snm(n::Int, m::Int, N::Int, ::Type{T}=Float64) where T = periodicleastsquarescoefficients(n, m, N, T)


IMPLEMENTED_FILTERS = CompactInfiniteVector[]
for n in 0:5
    local bn, bnm
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
