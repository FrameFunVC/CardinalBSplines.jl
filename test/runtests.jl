using CardinalBSplines, LinearAlgebra, Test, QuadGK
CREATE_README = false



function elementarypropsofsplinetest(T)
    T = real(T)
    tol = sqrt(eps(T))
    S = 20
    for N in 1:10
        f = BSpline(N-1,T)
        # Integral should be 1
        if !(T <: BigFloat)
            I,e = QuadGK.quadgk(f, 0, N, rtol = tol)
            @test I≈T(1)
        end
        # Infinite summation of shifted versions is 1
        xx = LinRange(T(N-1), T(N), S)[1:end-1]
        # xx = range(T(N-1), step=T(1)/T(S), length=S)
        g = zeros(T,length(xx))
        for k in 0:N-1
            g += map(x->f(x-k), xx)
        end
        @test g ≈ ones(T,length(g))  # (norm(g-1) < tol)
        # Two scale relation
        x = LinRange(T(-1), T(N+1), S)
        # x = range(T(-1), step=T(1)/T(S), length=S)
        g = zeros(T,length(x))
        for k in 0:N
            g += T(binomial(N,k))*map(x->f(2x-k), x)
        end
        g *= T(2)^(-N+1)
        G = map(f, x)
        @test g ≈ G
    end
end

function periodicbsplinetest(T)
    for N in 0:4
        period = T(N+1)
        x = LinRange(T(0),period,10)[1:end-1]
        PeriodicBSpline(N, period, T).(x) ≈ BSpline(N,T).(x)
        for k in -2:2
            @test PeriodicBSpline(N, period, T).(x) ≈ PeriodicBSpline(N,period,T).(x .+ k*period)
        end
        period = T(N+1)/T(3)
        x = LinRange(T(0),period,10)[1:end-1]
        for k in -2:2
            @test PeriodicBSpline(N, period, T).(x) ≈ PeriodicBSpline(N,period,T).(x .+ k*period)
        end
    end
end

function centeredbsplinestest(T)
    K = 20
    for N in 0:4
        x = LinRange(T(0)+eps(T), T(1), K)
        @test PeriodicCenteredBSpline(N, 10(N+1), T).(x) ≈ PeriodicCenteredBSpline(N, 10(N+1), T).(-x)
        @test CenteredBSpline(N, T).(x) ≈ CenteredBSpline(N, T).(-x)
    end
end

function test_spline_integration()
    T = Float64
    for N in 0:8
        @test squared_spline_integral(N) ≈ shifted_spline_integral(N,0)
        @test squared_spline_integral(N, Float64) ≈ shifted_spline_integral(N,0, Float64)
        @test squared_spline_integral(N) ≈ shifted_spline_integral(N,0, Float64)
        for t in 0:4
            f = BSpline(N, T)
            i = QuadGK.quadgk(x->f(x)*f(x+t),(-N-2:N+2)..., rtol=sqrt(eps(T)))[1]
            @test abs(i-shifted_spline_integral(N,t))<10*sqrt(eps(T))
            @test abs(i-shifted_spline_integral(N,t,Float64))<10*sqrt(eps(T))
            @test abs(i-shifted_spline_integral(N,t,BigFloat))<10*sqrt(eps(T))
        end
    end
end


function derivative_test(T)
    t = LinRange(-10,10,100); t2 = (t[1:end-1]+t[2:end])/2
    for i in 2:10
        y1 = BSpline(i,T).(t)
        y2 = BSplineDiff(i,T).(t2)
        @test norm(diff(y1)./diff(t).-y2) < .1
    end
    for i in 2:10
        y1 = PeriodicBSpline(i, 5, T).(t)
        y2 = PeriodicBSplineDiff(i, 5, T).(t2)
        @test norm(diff(y1)./diff(t).-y2) < .1
    end

    for t in [LinRange(0,1,100), LinRange(-1,0,100), LinRange(1,2,100)]
        t2 = (t[1:end-1]+t[2:end])/2
        y1 = BSpline(1,T).(t)
        y2 = BSplineDiff(1,T).(t2)
        @test norm(diff(y1)./diff(t).-y2) < .1
    end

    t0 = LinRange(-1,7,100)
    for D in 0:6
        @test norm(BSplineDiff(D,D+1,T).(t0))==0
    end

    t0 = -2.25:1/(1<<4):10
    t1 = (t0[1:end-1]+t0[2:end])/2
    t2 = (t1[1:end-1]+t1[2:end])/2
    t3 = (t2[1:end-1]+t2[2:end])/2
    t4 = (t3[1:end-1]+t3[2:end])/2
    t5 = (t4[1:end-1]+t4[2:end])/2

    for D in 3:6
        y0 = BSplineDiff(D,0,T).(t0)
        y2 = BSplineDiff(D,2,T).(t2)
        @test .5 > norm(diff(diff(y0)./diff(t0))./diff(t1) - y2)
    end

    for D in 4:6
        y0 = BSplineDiff(D,0,T).(t0)
        y2 = BSplineDiff(D,3,T).(t3)
        @test .5 > norm(  diff(diff(diff(y0)./diff(t0))./diff(t1))./diff(t2) - y2)
    end

    for D in 5:6
        y0 = BSplineDiff(D,0,T).(t0)
        y2 = BSplineDiff(D,4,T).(t4)
        @test .5 > norm(  diff(diff(diff(diff(y0)./diff(t0))./diff(t1))./diff(t2))./diff(t3) - y2)
    end

    for D in 6:8
        y0 = BSplineDiff(D,0,T).(t0)
        y2 = BSplineDiff(D,5,T).(t5)
        @test .5 > norm(  diff(diff(diff(diff(diff(y0)./diff(t0))./diff(t1))./diff(t2))./diff(t3))./diff(t4) - y2)
    end

    for D in 1:6
        @test Int.(BSplineDiff(D, D, T).(-.5:D+2)) == vcat(0,[(-1)^k*binomial(D,k) for k in 0:D],0)
    end
end


function allocation_test()
    for T in (BSpline, PeriodicBSpline, CenteredBSpline, PeriodicCenteredBSpline)
        for K in 0:10
            f = T(K)
            for x in LinRange(-10,10,30)
                f(x)
            end
            for x in LinRange(-10,10,30)
                @test @timed(f(x))[3] <= 32
            end
        end
    end
    for T in (BSplineDiff,)
        for K in 0:10, D in 1:K
            f = T(K,D)
            for x in LinRange(-10,10,2)
                f(x)
            end
            for x in LinRange(-10,10,30)
                @test @timed(f(x))[3] <= 32
            end
        end
    end
end

P = 80
for T in [Float64, BigFloat]
    @testset "$(rpad("Elementary properties",P))" begin
      elementarypropsofsplinetest(T)
    end
    @testset "$(rpad("periodic B splines",P))"  begin
      periodicbsplinetest(T)
    end
    @testset "$(rpad("centered B splines",P))"  begin
      centeredbsplinestest(T)
    end
    @testset "$(rpad("derivative of B splines",P))"  begin
      derivative_test(T)
    end
end
@testset "$(rpad("integration of B splines",P))"  begin
    test_spline_integration()
end
@testset "$(rpad("allocation of evaluation of B splines",P))"  begin
    allocation_test()
end

include("filter_test.jl")

if CREATE_README
    cd(normpath((splitdir(@__FILE__))[1]*"/.."))
    try
        println("Create README.md")
        run(`jupyter nbconvert --execute --to markdown --output README.md notebooks/README.ipynb`)
        run(`mv notebooks/README.md .`)
        run(`mv notebooks/README_files/ .`)
    catch
        nothing
    end
end
println("All tests succeeded")
