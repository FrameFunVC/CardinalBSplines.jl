export GaussSpline
struct GaussSpline{K,S} end

(GaussSpline)(K::Int, S=Float64) = (GaussSpline){K}(S)
(GaussSpline){K}(::Type{S}) where {K,S} = (GaussSpline){K,S}()
(::(GaussSpline){K,S})(x) where {K,S} = evaluate_centered_gauss_BSpline(Val{K}(), convert(S,x), S)

function evaluate_centered_gauss_BSpline(::Val{K}, x, ::Type{T}) where {K,T}
    sqrt(T(6)/T(pi)/T(K+1))*exp(-T(6)*x^2/T(K+1))
end

abstract type AbstractBSplineDerivative{K,D,S} end

abstract type AbstractBSpline{K,S} <: AbstractBSplineDerivative{K,0,S} end




for (T, f) in ((:BSpline, :evaluate_BSpline),
                (:CenteredBSpline, :evaluate_centered_BSpline))
    @eval begin
        struct ($T){K,S} <:AbstractBSpline{K,S}
            function ($T){K,S}() where {K,S}
                # @assert K isa Integer && K >= 0
                # @assert S <: Real
                new{K,S}()
            end
        end

        ($T)() = ($T){3}()
        ($T)(::Type{S}) where {S} = ($T){3}(S)
        ($T)(::Val{K}, args...) where {K} = ($T){K}(args...)
        ($T){K}(::Type{S}=Float64) where {K,S} = ($T){K,S}()
        (::($T){K,S})(x) where {K,S} = $f(Val{K}(), convert(S,x), S)
        # Warning: the one below is not type stable, but it is convenient
        ($T)(K::Int, ::Type{S}=Float64) where {S} = ($T){K,S}()
    end
end

export BSpline, CenteredBSpline

@doc """
    BSpline(m, T=Float64)

Construct a B-spline of order `m` and type `T`


The support of the B-spline is [0,m+1].
"""
BSpline(m, T)

@doc """
    CenteredBSpline(m, T=Float64)

Construct a centered B-spline of order `m` and type `T`


The support of the centered B-spline is s \$\\left[-\\tfrac{m+1}{2},\\tfrac{m+1}{2}\\right]\$.
"""
CenteredBSpline(m, T)


for (T, f) in ((:PeriodicBSpline, :evaluate_periodic_BSpline),
                (:PeriodicCenteredBSpline, :evaluate_periodic_centered_BSpline))
    @eval begin
        struct ($T){K,S} <:AbstractBSpline{K,S}
            period  :: S
            function ($T){K,S}(period) where {K,S}
                # @assert K isa Integer && K >= 0
                # @assert S <: Real
                # @assert period > 0
                new{K,S}(convert(S, period))
            end
        end

        ($T)() = ($T){3}()
        ($T)(::Val{K}, args...) where {K} = ($T){K}(args...)
        ($T)(::Type{S}) where {S} = ($T){3}(S)
        ($T){K}(::Type{S}) where {K,S} = ($T){K}(1, S)
        ($T){K}(period = 1, ::Type{S}=Float64) where {K,S} = ($T){K,S}(period)
        (spline::($T){K,S})(x) where {K,S} = $f(Val{K}(), convert(S,x), spline.period, S)
        # Warning: the one below is not type stable, but it is convenient
        ($T)(K::Int, period=1, S=Float64) = ($T){K}(period, S)
    end
end

export PeriodicBSpline, PeriodicCenteredBSpline

@doc """
    PeriodicBSpline(m, p=1, S=Float64)

Construct a periodic B-spline of order `m`, period `p`, and type `S`

See also `BSpline`
"""
PeriodicBSpline(m, p, S)

@doc """
    PeriodicCenteredBSpline(m, p=1, S=Float64)

Construct a periodic centered B-spline of order `m` and type `S`.

See also `CenteredBSpline`
"""
PeriodicCenteredBSpline(m, p, S)

for (T, f) in ((:BSplineDiff, :evaluate_BSpline_derivative),
                (:CenteredBSplineDiff, :evaluate_centered_BSpline_derivative))
    @eval begin
        struct ($T){K,D,S}<:AbstractBSplineDerivative{K,D,S}
            function ($T){K,D,S}() where {K,D,S}
                # @assert K isa Integer && K >= 0
                # @assert D isa Integer && D >= 0
                # @assert S <: Real
                new{K,D,S}()
            end
        end

        # First, parse spline degree K
        ($T)() = ($T){3}()
        ($T)(::Val{K}, args...) where {K} = ($T){K}(args...)
        # Next, parse the differentiation order
        # - the order is omitted, only S is given: use D=1
        ($T){K}(::Type{S}) where {K,S} = ($T){K}(Val{1}(), S)
        # - order is specified as Val{D} or nothing is given
        ($T){K}(d::Val{D}=Val{1}(), args...) where {K,D} = ($T){K,D}(args...)
        # Finally, parse S
        ($T){K,D}(::Type{S}=Float64) where {K,D,S} = ($T){K,D,S}()
        (::($T){K,D,S})(x) where {K,D,S} = $f(Val{K}(), Val{D}(), convert(S,x), S)
        # Warning: the ones below are not type stable, but they are convenient
        ($T)(K::Int, args...) = ($T){K}(args...)
        ($T){K}(D::Int, args...) where {K} = ($T){K,D}(args...)
    end
end


export BSplineDiff, CenteredBSplineDiff

@doc """
    BSplineDiff(m, d, S=Float64)

Construct the `d`th derivative of a B-spline of order `m` and type `S`

See also `BSpline`
"""
BSplineDiff(m, d, S)

@doc """
    CenteredBSplineDiff(m, d, S=Float64)

Construct the `d`th derivative of a centered B-spline of order `m` and type `S`.

See also `CenteredBSpline`
"""
CenteredBSplineDiff(m, d, S)



for (T, f) in ((:PeriodicBSplineDiff, :evaluate_periodic_BSpline_derivative),
                (:PeriodicCenteredBSplineDiff, :evaluate_periodic_centered_BSpline_derivative))
    @eval begin
        struct ($T){K,D,S} <:AbstractBSplineDerivative{K,D,S}
            period  :: S
            function ($T){K,D,S}(period) where {K,D,S}
                # @assert K isa Integer && K >= 0
                # @assert D isa Integer && D >= 0
                # @assert S <: Real
                # @assert period > 0
                new{K,D,S}(convert(S, period))
            end
        end

        # The order of the arguments is: K, period, D, S.

        # First: parse the spline degree and invoke $T{K}
        ($T)() = ($T){3}()
        ($T)(::Val{K}, args...) where {K} = ($T){K}(args...)
        # Next: parse the period
        # - the period is omitted (obvious because Val{D} is given): add default
        ($T){K}(d::Val{D}, args...) where {K,D} = ($T){K}(1, d, args...)
        # - assume defaults and invoke $T{K,D}
        ($T){K}(period::Number=1, d::Val{D}=Val{1}(), args...) where {K,D} = ($T){K,D}(period, args...)
        # - the case where period and S are given, but D is not
        ($T){K}(period::Number, ::Type{S}) where {K,S} = ($T){K}(period, Val{1}(), S)
        # - only S is given
        ($T){K}(::Type{S}) where {K,S} = ($T){K}(1, Val{1}(), S)
        # Finally, invoke $T{K,D,S}
        ($T){K,D}(::Type{S}=Float64) where {K,D,S} = ($T){K,D}(1, S)
        ($T){K,D}(period::Number, ::Type{S}=Float64) where {K,D,S} = ($T){K,D,S}(period)
        (spline::($T){K,D,S})(x) where {K,S,D} = $f(Val{K}(), Val{D}(), convert(S,x), spline.period, S)
        # Warning: the ones below are not type stable, but they are convenient
        ($T)(K::Int, args...) = ($T){K}(args...)
        ($T){K}(period::Number, D::Int, args...) where {K} = ($T){K}(period, Val{D}(), args...)
    end
end

export PeriodicBSplineDiff, PeriodicCenteredBSplineDiff
@doc """
    PeriodicBSplineDiff(m, p, d, S=Float64)

Construct the `d`th derivative of a periodic B-spline of order `m` , period `p` and type `S`

See also `PeriodicBSpline`
"""
PeriodicBSplineDiff(m, p, d, S)

@doc """
    PeriodicCenteredBSplineDiff(m, p, d, S=Float64)

Construct the `d`th derivative of a periodic, centered B-spline of order `m`, period `p`, and type `S`.

See also `PeriodicCenteredBSpline`
"""
PeriodicCenteredBSplineDiff(m, p, d, S)
