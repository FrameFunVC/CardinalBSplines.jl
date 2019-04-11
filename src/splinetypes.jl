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
