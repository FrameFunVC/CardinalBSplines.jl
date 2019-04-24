
abstract type AbstractBSplineDerivative{K,D,S} end

abstract type AbstractBSpline{K,S} <: AbstractBSplineDerivative{K,0,S} end




for (T, f) in ((:BSpline, :evaluate_Bspline),
                (:CenteredBSpline, :evaluate_centered_BSpline))
    @eval begin
        struct ($T){K,S} <:AbstractBSpline{K,S}
            function ($T){K,S}() where {K,S}
                # @assert K isa Integer && K >= 0
                # @assert S <: Real
                new{K,S}()
            end
        end

        ($T)(K::Int, S=Float64) = ($T){K}(S)
        ($T){K}(::Type{S}) where {K,S} = ($T){K,S}()
        (::($T){K,S})(x) where {K,S} = $f(Val{K}(), convert(S,x), S)
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


for (T, f) in ((:PeriodicBSpline, :evaluate_periodic_Bspline),
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

        ($T)(K::Int, period=1, S=Float64) = ($T){K}(period, S)
        ($T){K}(period, ::Type{S}) where {K,S} = ($T){K,S}(period)
        (spline::($T){K,S})(x) where {K,S} = $f(Val{K}(), convert(S,x), spline.period, S)
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

for (T, f) in ((:BSplineDiff, :evaluate_Bspline_derivative),
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

        ($T)(K::Int, S=Float64) = ($T)(K, 1, S)
        ($T)(K::Int, D::Int, S=Float64) = ($T){K,D}(S)
        ($T){K,D}(::Type{S}) where {K,D,S} = ($T){K,D,S}()
        (::($T){K,D,S})(x) where {K,D,S} = $f(Val{K}(), Val{D}(), convert(S,x), S)
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



for (T, f) in ((:PeriodicBSplineDiff, :evaluate_periodic_Bspline_derivative),
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
        ($T)(K::Int, period=1, S=Float64) = ($T)(K, period, 1, S)
        ($T)(K::Int, period, D::Int, S=Float64) = ($T){K,D}(period, S)
        ($T){K,D}(period, ::Type{S}) where {K,D,S} = ($T){K,D,S}(period)
        (spline::($T){K,D,S})(x) where {K,S,D} = $f(Val{K}(), Val{D}(), convert(S,x), spline.period, S)
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
