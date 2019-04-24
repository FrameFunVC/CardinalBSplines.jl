export squared_spline_integral, shifted_spline_integral
"""
    squared_spline_integral(N)

Calculate ∫\$[B_N(x)]^2\$ dx, with \$B_N\$ the B-spline of degree N
"""
squared_spline_integral(N::Int, ::Type{T}=BigInt) where T<:Real =
    squared_spline_integral(N, zero(T))

"""
    shifted_spline_integral(m, t)

Calculate ∫ \$B_m(x) B_m(x-t)\$ dx, with \$B_m\$ the B-spline of degree m
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
