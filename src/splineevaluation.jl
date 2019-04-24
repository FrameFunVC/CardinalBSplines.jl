
# Implementation of cardinal B splines of degree N
"""
    evaluate_periodic_centered_BSpline(::Val{N}, x, period, ::Type{T})

Evaluate the periodic centered B-spline of order `N`, type `T` and period `p` in `x`.
"""
evaluate_periodic_centered_BSpline(::Val{N}, x::Real, period::Real, ::Type{T}) where {T<:Real,N} =
    evaluate_periodic_Bspline(Val{N}(), x+T(N+1)/2, period, T)

"""
    evaluate_centered_BSpline(::Val{N}, x, period, ::Type{T})

Evaluate the centered B-spline of order `N` and type `T` in `x`.
"""
evaluate_centered_BSpline(::Val{N}, x::Real, ::Type{T}) where {T<:Real,N} =
    evaluate_Bspline(Val{N}(), x+T(N+1)/2, T)

function periodize(x::Real, period::Real)
    x -= period*fld(x, period)
    x >= period && (x -= period)
    # @assert(0<= x < period)
    x
end

"""
    evaluate_periodic_BSpline(::Val{N}, x, period, ::Type{T})

Evaluate the periodic B-spline of order `N`, type `T` and period `p` in `x`.
"""
function evaluate_periodic_Bspline(::Val{N}, x::Real, period::Real, ::Type{T}) where {N,T<:Real}
    # @assert period > 0
    x = periodize(x, period)
    res = T(0)
    for k in 0:floor(Int, (N+1-x)/period)
        res += evaluate_Bspline(Val{N}(), x+period*k, T)
    end
    res
end

"""
    evaluate_BSpline(::Val{N}, x, ::Type{T})

Evaluate the periodic centered B-spline of order `N`, and type `T` in `x`.
"""
function evaluate_Bspline(::Val{N}, x::Real, ::Type{T}) where {N,T<:Real}
    (T(x)/T(N)*evaluate_Bspline(Val{N-1}(), x, T) +
        (T(N+1)-T(x))/T(N)*evaluate_Bspline(Val{N-1}(), x-1, T))
end

evaluate_Bspline(::Val{0}, x::Real, ::Type{T}) where {T<:Real} = (0 <= x < 1) ? T(1) : T(0)

function evaluate_Bspline(::Val{1}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return T(x)
    elseif (1 <= x < 2)
        return T(2) - T(x)
    else
        return T(0)
    end
end

@eval function evaluate_Bspline(::Val{2}, x::Real, ::Type{T}) where {T<:Real}
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

@eval function evaluate_Bspline(::Val{3}, x::Real, ::Type{T}) where {T<:Real}
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

@eval function evaluate_Bspline(::Val{4}, x::Real, ::Type{T}) where {T<:Real}
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

"""
    evaluate_BSpline_derivative(::Val{N}, ::Var{D}, x, ::Type{T})

Evaluate the `D`th derivative of the B-spline of order `N`, and type `T` in `x`.
"""
evaluate_Bspline_derivative(::Val{N}, ::Val{0}, x::Real, ::Type{T}) where {N,T<:Real} =
    evaluate_Bspline(Val{N}(), x, T)

evaluate_Bspline_derivative(::Val{0}, ::Val{1}, x::Real, ::Type{T}) where {T<:Real} =
    T(0)

evaluate_Bspline_derivative(::Val{0}, ::Val{K}, x::Real, ::Type{T}) where {K,T<:Real} =
    T(0)
evaluate_Bspline_derivative(::Val{0}, ::Val{0}, x::Real, ::Type{T}) where {K,T<:Real} =
    evaluate_Bspline(Val{0}(), x, T)


function evaluate_Bspline_derivative(::Val{1}, ::Val{1}, x::Real, ::Type{T}) where {T<:Real}
    if (0 <= x < 1)
        return T(1)
    elseif (1 <= x < 2)
        return - T(1)
    else
        return T(0)
    end
end

@eval function evaluate_Bspline_derivative(::Val{2}, ::Val{1}, x::Real, ::Type{T}) where {T<:Real}
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

@eval function evaluate_Bspline_derivative(::Val{3}, ::Val{1}, x::Real, ::Type{T}) where {T<:Real}
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

@eval function evaluate_Bspline_derivative(::Val{4}, ::Val{1}, x::Real, ::Type{T}) where {T<:Real}
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


function evaluate_Bspline_derivative(::Val{1}, ::Val{2}, x::Real, ::Type{T}) where {T<:Real}
    T(0) # +1 at 0, -2 at 1, +1 at 2
end

function evaluate_Bspline_derivative(::Val{2}, ::Val{2}, x::Real, ::Type{T}) where {T<:Real}
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

@eval function evaluate_Bspline_derivative(::Val{3}, ::Val{2}, x::Real, ::Type{T}) where {T<:Real}
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

@eval function evaluate_Bspline_derivative(::Val{4}, ::Val{2}, x::Real, ::Type{T}) where {T<:Real}
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

function evaluate_Bspline_derivative(::Val{2}, ::Val{3}, x::Real, ::Type{T}) where {T<:Real}
    T(0)# jump 1 at 0, -3 at 1, +3, at 2 -1 at 3
end

function evaluate_Bspline_derivative(::Val{3}, ::Val{3}, x::Real, ::Type{T}) where {T<:Real}
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

@eval function evaluate_Bspline_derivative(::Val{4}, ::Val{3}, x::Real, ::Type{T}) where {T<:Real}
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

function evaluate_Bspline_derivative(::Val{3}, ::Val{4}, x::Real, ::Type{T}) where {T<:Real}
    T(0) # jump 1 at 0, -4 at 2, +6 at 3, -4 at 4, +1, at 5
end

function evaluate_Bspline_derivative(::Val{4}, ::Val{4}, x::Real, ::Type{T}) where {T<:Real}
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

function evaluate_Bspline_derivative(::Val{4}, ::Val{5}, x::Real, ::Type{T}) where {T<:Real}
    T(0) # jump 1 at 0, -5 at 1, +10 at 2 -10 at 3, +5 at 4, -1 at 5.
end

function evaluate_Bspline_derivative(::Val{N}, ::Val{N}, x::Real, ::Type{T}) where {N,T<:Real}
    xfl = floor(Int, x)
    if xfl < 0 || xfl > N
        T(0)
    else
        T((-1)^xfl*binomial(N,xfl))
    end
end

# reduce degree
function evaluate_Bspline_derivative(::Val{N}, ::Val{K}, x::Real, ::Type{T}) where {N,K,T}
    T(K)/T(N)*(evaluate_Bspline_derivative(Val{N-1}(), Val{K-1}(), x, T) - evaluate_Bspline_derivative(Val{N-1}(), Val{K-1}(), x-1, T)) +
        (T(N+1)-T(x))/T(N)*evaluate_Bspline_derivative(Val{N-1}(), Val{K}(), x-1, T) + T(x)/T(N)*evaluate_Bspline_derivative(Val{N-1}(), Val{K}(), x, T)
end

"""
    evaluate_periodic_BSpline_derivative(::Val{N}, ::Var{D}, x, period, ::Type{T})

Evaluate the `D`th derivative of the periodic B-spline of order `N`, type `T` and period `p` in `x`.
"""
function evaluate_periodic_Bspline_derivative(::Val{N}, ::Val{K}, x::Real, period::Real, ::Type{T}) where {N,K,T<:Real}
    x = periodize(x, period)
    res = T(0)
    for k in 0:floor(Int, (N+1-x)/period)
        res += evaluate_Bspline_derivative(Val{N}(), Val{K}(), x+period*k, T)
    end
    res
end

"""
    evaluate_periodic_centered_BSpline_derivative(::Val{N}, ::Var{D}, x, period, ::Type{T})

Evaluate the `D`th derivative of the periodic centered B-spline of order `N`, type `T` and period `p` in `x`.
"""
evaluate_periodic_centered_BSpline_derivative(::Val{N}, ::Val{K}, x::Real, period::Real, ::Type{T}) where {N,K,T<:Real} =
    evaluate_periodic_BSpline_derivative(Val{N}(), Val{K}(), x+T(N+1)/2, period, T)


"""
    evaluate_centered_BSpline_derivative(::Val{N}, ::Var{D}, x, ::Type{T})

Evaluate the `D`th derivative of the centered B-spline of order `N`, type `T` and period `p` in `x`.
"""
evaluate_centered_BSpline_derivative(::Val{N}, ::Val{K}, x::Real, ::Type{T}) where {N,K,T<:Real} =
    evaluate_BSpline_derivative(Val{N}(), Val{K}(), x+T(N+1)/2, T)
