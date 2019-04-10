# CardinalBSplines.jl
module CardinalBSplines
using SpecialFunctions: gamma

export evaluate_Bspline, evaluate_periodic_Bspline, Degreek, evaluate_symmetric_periodic_Bspline, squared_spline_integral, shifted_spline_integral
export diff_evaluate_Bspline, diff_evaluate_periodic_Bspline
"""
  Calculates ∫[B_N(x)]^2 dx, with B_N the cardinal bspline of degree N
"""
squared_spline_integral(N::Int, ::Type{T}=BigInt) where T<:Real =
    squared_spline_integral(N, zero(T))

"""
  Calculates ∫ B_m(x) B_m(x-t) dx, with B_m the cardinal bspline of degree m
"""
shifted_spline_integral(m::Int,t::Int, ::Type{T}=BigInt) where T<:Real =
    shifted_spline_integral(m, t, zero(T))

function squared_spline_integral(N::Int, ::T) where T<:AbstractFloat
    m = convert(T,N+1)
    S = zero(T)
    for k in zero(T):m-1
        S += binomial(2m,k)*(-1)^k*(2*(m-k))^(2m-1)
    end
    S/(convert(T,2)^(2m-1)*factorial(2m-1))
end

Base.binomial(n::T,k::T) where T<:Real = gamma(n+1)/gamma(k+1)/gamma(n-k+1)
Base.factorial(n::T) where T<:Real = gamma(n+1)

function shifted_spline_integral(m::Int,t::Int, ::T) where T<:AbstractFloat
    @assert t >= 0
    E = convert(T,t)
    M = convert(T, m)
    I = zero(T)
    O = one(T)
    if t > m
        return zero(T)
    end
    for k in zero(T):M+O
        for j in zero(T):M+O
            S = zero(T)
            if j <= k-E
                for i in zero(T):M
                    S += binomial(M,i)*(k-j-E)^(m-i)*(M+O-k)^(i+M+O)/(i+m+O)
                end
                I += (-O)^(j+k)*binomial(M+O,k)*binomial(M+O,j)*S
            else
                for i in zero(T):M
                    S += binomial(M,i)*(E+j-k)^(M-i)*(M+O-j-E)^(i+M+O)/(i+M+O)
                end
                I += (-O)^(j+k)*binomial(M+O,k)*binomial(M+O,j)*S
            end
        end
    end
    convert(T,I/(factorial(M)^convert(T,2)))
end

"""
  Calculates ∫[B_N(x)]^2 dx, with B_N the cardinal bspline of degree N
"""
function squared_spline_integral(N::Int, ::T) where T<:Integer
    m = convert(T, N+1)
    if m > 7 && sizeof(T) <= 8
        error("You need a higher precision, try `squared_spline_integral($N, BigInt)`")
    end
    S = zero(Rational{T})
    for k in 0:m-1
        S += binomial(2m,k)*(-1)^k*(2*(m-k))^(2m-1)
    end
    S//(2^(2m-1)*factorial(2m-1))
end

"""
  Calculates ∫ B_m(x) B_m(x-t) dx, with B_m the cardinal bspline of degree m
"""
function shifted_spline_integral(m::Int,t::Int, ::T) where T<:Integer
    # @assert t >= 0
    t = abs(t)
    if t == 0
        return squared_spline_integral(m, zero(T))
    end
    if t > m
        return convert(Rational{T},0)
    end
    m = convert(T,m)
    I = convert(Rational{T},0)
    if m > 7 && sizeof(T) <= 8
        error("You need a higher precision, try `shifted_spline_integral($m, $t, BigInt)`")
    end
    for k in 0:m+1
        for j in 0:m+1
            if j <= k-t
                S = zero(Rational{T})
                for i in 0:m
                    S += binomial(m,i)*(k-j-t)^(m-i)*(m+1-k)^(i+m+1)//(i+m+1)
                end
                I = I + convert(Rational{T},-1)^(j+k)*binomial(m+1,k)*binomial(m+1,j)*S
            else
                S = zero(Rational{T})
                for i in 0:m
                    S += binomial(m,i)*(t+j-k)^(m-i)*(m+1-j-t)^(i+m+1)//(i+m+1)
                end
                I = I + convert(Rational{T},-1)^(j+k)*binomial(m+1,k)*binomial(m+1,j)*S
            end
        end
    end
    I//(factorial(m)^2)
end

# Implementation of cardinal B splines of degree N
Degree{N} = Val{N}
function evaluate_symmetric_periodic_Bspline(N::Int, x::Real, period::Real, ::Type{T}) where {T<:Real}
    evaluate_periodic_Bspline(N, x+T(N+1)/2, period, T)
end

function periodize(x::Real, period::Real)
    x -= period*fld(x, period)
    x >= period && (x -= period)
    @assert(0<= x < period)
    x
end

evaluate_periodic_Bspline(N::Int, x::Real, period::Real, ::Type{T}) where {T} =
    evaluate_periodic_Bspline(Val{N}, x, period, T)

function evaluate_periodic_Bspline(::Type{Degree{N}}, x::Real, period::Real, ::Type{T}) where {N,T<:Real}
    @assert period > 0
    x = periodize(x, period)
    res = T(0)
    for k in 0:floor(Int, (N+1-x)/period)
        res += evaluate_Bspline(Val{N}, x+period*k, T)
    end
    res
end

evaluate_Bspline(N::Int, x::Real, ::Type{T}) where {T<:Real} = evaluate_Bspline(Degree{N}, x, T)

evaluate_Bspline(::Type{Degree{N}}, x::Real, ::Type{T}) where {N,T<:Real} =
    _evaluate_Bspline(N, x, T)
function _evaluate_Bspline(N::Int, x::Real, ::Type{T}) where {T<:Real}
    if N == 5
        (T(x)/T(N)*evaluate_Bspline(Val{4}, x, T) +
            (T(N+1)-T(x))/T(N)*evaluate_Bspline(Val{4}, x-1, T))
    else
        (T(x)/T(N)*_evaluate_Bspline(N-1, x, T) +
            (T(N+1)-T(x))/T(N)*_evaluate_Bspline(N-1, x-1, T))
    end
end

evaluate_Bspline(::Type{Degree{0}}, x::Real, ::Type{T}) where {T<:Real} = (0 <= x < 1) ? T(1) : T(0)

function evaluate_Bspline(::Type{Degree{1}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return T(x)
    elseif (1 <= x < 2)
        return T(2) - T(x)
    else
        return T(0)
    end
end

@eval function evaluate_Bspline(::Type{Degree{2}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return @evalpoly(T(x),T(0), T(0), T(1/2))
    elseif (1 <= x < 2)
        return @evalpoly(T(x),T(-3/2), T(3), T(-1))
    elseif (2 <= x < 3)
        return @evalpoly(T(x),T(9//2), T(-3), T(1//2))
    else
        return T(0)
    end
end

@eval function evaluate_Bspline(::Type{Degree{3}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(0), T(0), T(1//6))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(2//3), T(-2), T(2), T(-1//2))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(-22//3), T(10), T(-4), T(1//2))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(32//3), T(-8), T(2), T(-1//6))
    else
        return T(0)
    end
end

@eval function evaluate_Bspline(::Type{Degree{4}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(0), T(0), T(0), T(1//24))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(-5//24), T(5//6), T(-5/4), T(5//6), T(-1//6))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(155//24), T(-25//2), T(35//4), T(-5//2), T(1//4))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(-655//24), T(65//2), T(-55//4), T(5//2), T(-1//6))
    elseif (4 <= x < 5)
        return @evalpoly(T(x), T(625//24), -T(125//6), T(25//4), T(-5//6), T(1//24))
    else
        return T(0)
    end
end

diff_evaluate_Bspline(::Type{Degree{N}}, ::Type{Val{0}}, x::Real, ::Type{T}) where {N,T<:Real} =
    evaluate_Bspline(N, x, T)

function diff_evaluate_Bspline(::Type{Val{0}}, ::Type{Val{1}}, x::Real, ::Type{T}) where {T<:Real}
    # if (x == 0 || x == 1)
    #     return NaN
    # else
    #     return T(0)
    # end
    T(0)
end

function diff_evaluate_Bspline(::Type{Degree{1}}, ::Type{Val{1}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return T(1)
    elseif (1 <= x < 2)
        return - T(1)
    else
        return T(0)
    end
end

@eval function diff_evaluate_Bspline(::Type{Degree{2}}, ::Type{Val{1}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return @evalpoly(T(x),T(0), T(1))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(3), T(-2))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(-3), T(1))
    else
        return T(0)
    end
end

@eval function diff_evaluate_Bspline(::Type{Degree{3}}, ::Type{Val{1}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(0), T(1//2))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(-2), T(4), T(-3//2))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(10), T(-8), T(3//2))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(-8), T(4), T(-1//2))
    else
        return T(0)
    end
end

@eval function diff_evaluate_Bspline(::Type{Degree{4}}, ::Type{Val{1}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(0), T(0), T(1//6))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(5//6), T(-5/2), T(5//2), T(-2//3))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(-25//2), T(35//2), T(-15//2), T(1))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(65//2), T(-55//2), T(15//2), T(-2//3))
    elseif (4 <= x < 5)
        return @evalpoly(T(x), -T(125//6), T(25//2), T(-5//2), T(1//6))
    else
        return T(0)
    end
end


function diff_evaluate_Bspline(::Type{Degree{1}}, ::Type{Val{2}}, x::Real, ::Type{T}) where {T<:Real}
    T(0) # +1 at 0, -2 at 1, +1 at 2
end

function diff_evaluate_Bspline(::Type{Degree{2}}, ::Type{Val{2}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return T(1)
    elseif (1 <= x < 2)
        return T(-2)
    elseif (2 <= x < 3)
        return T(1)
    else
        return T(0)
    end
end

@eval function diff_evaluate_Bspline(::Type{Degree{3}}, ::Type{Val{2}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(1))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(4), T(-3))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(-8), T(3))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(4), T(-1))
    else
        return T(0)
    end
end

@eval function diff_evaluate_Bspline(::Type{Degree{4}}, ::Type{Val{2}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(0), T(1//2))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(-5/2), T(5), T(-2))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(35//2), T(-15), T(3))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(-55//2), T(15), T(-2))
    elseif (4 <= x < 5)
        return @evalpoly(T(x), T(25//2), T(-5), T(1//2))
    else
        return T(0)
    end
end

function diff_evaluate_Bspline(::Type{Degree{2}}, ::Type{Val{3}}, x::Real, ::Type{T}) where {T<:Real}
    T(0)# jump 1 at 0, -3 at 1, +3, at 2 -1 at 3
end

function diff_evaluate_Bspline(::Type{Degree{3}}, ::Type{Val{3}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return T(1)
    elseif (1 <= x < 2)
        return T(-3)
    elseif (2 <= x < 3)
        return  T(3)
    elseif (3 <= x < 4)
        return T(-1)
    else
        return T(0)
    end
end

@eval function diff_evaluate_Bspline(::Type{Degree{4}}, ::Type{Val{3}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return @evalpoly(T(x), T(0), T(1))
    elseif (1 <= x < 2)
        return @evalpoly(T(x), T(5), T(-4))
    elseif (2 <= x < 3)
        return @evalpoly(T(x), T(-15), T(6))
    elseif (3 <= x < 4)
        return @evalpoly(T(x), T(15), T(-4))
    elseif (4 <= x < 5)
        return @evalpoly(T(x), T(-5), T(1))
    else
        return T(0)
    end
end

function diff_evaluate_Bspline(::Type{Degree{3}}, ::Type{Val{4}}, x::Real, ::Type{T}) where {T<:Real}
    T(0) # jump 1 at 0, -4 at 2, +6 at 3, -4 at 4, +1, at 5
end

function diff_evaluate_Bspline(::Type{Degree{4}}, ::Type{Val{4}}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return T(1)
    elseif (1 <= x < 2)
        return T(-4)
    elseif (2 <= x < 3)
        return T(6)
    elseif (3 <= x < 4)
        return T(-4)
    elseif (4 <= x < 5)
        return  T(1)
    else
        return T(0)
    end
end

function diff_evaluate_Bspline(::Type{Degree{4}}, ::Type{Val{5}}, x::Real, ::Type{T}) where {T<:Real}
    T(0) # jump 1 at 0, -5 at 1, +10 at 2 -10 at 3, +5 at 4, -1 at 5.
end

function diff_evaluate_Bspline(::Type{Val{N}}, ::Type{Val{N}}, x::Real, ::Type{T}) where {N,T<:Real}
    xfl = floor(Int, x)
    if xfl < 0 || xfl > N
        T(0)
    else
        T((-1)^xfl*binomial(N,xfl))
    end
end

# reduce degree
function diff_evaluate_Bspline(::Type{Degree{N}}, ::Type{Val{K}}, x::Real, T::Type) where {N,K}
  T(K)/T(N)*(diff_evaluate_Bspline(Degree{N-1}, Val{K-1}, x, T) - diff_evaluate_Bspline(Degree{N-1}, Val{K-1}, x-1, T)) +
      (T(N+1)-T(x))/T(N)*diff_evaluate_Bspline(Degree{N-1}, Val{K}, x-1, T) + T(x)/T(N)*diff_evaluate_Bspline(Degree{N-1}, Val{K}, x, T)
end

diff_evaluate_Bspline(::Type{Degree{0}}, ::Type{Val{K}}, x::Real, T::Type) where {K} = T(0)

diff_evaluate_Bspline(N::Int, x::Real, ::Type{T}) where {T<:Real} = diff_evaluate_Bspline(Degree{N}, Val{1}, x, T)

diff_evaluate_periodic_Bspline(N::Int, x::Real, period::Real, ::Type{T}) where {T<:Real} =
    diff_evaluate_periodic_Bspline(Val{N}, Val{1}, x, period, T)

function diff_evaluate_periodic_Bspline(::Type{Val{N}}, ::Type{Val{K}}, x::Real, period::Real, ::Type{T}) where {N,K,T<:Real}
    @assert period > 0
    x = periodize(x, period)
    res = T(0)
    for k in 0:floor(Int, (N+1-x)/period)
        res += diff_evaluate_Bspline(Val{N}, Val{K}, x+period*k, T)
    end
    res
end

end # module Cardinal_b_splines
