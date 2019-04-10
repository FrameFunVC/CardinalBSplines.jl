using CardinalBSplines
using QuadGK

CREATE_README = false
if VERSION < v"0.7-"
    using Base.Test
    CREATE_README = true
else
    using LinearAlgebra
    using Test
    linspace(a,b,c) = range(a, stop=b, length=c)
end



function elementarypropsofsplinetest(T)
    T = real(T)
    tol = sqrt(eps(T))
    S = 20
    for N in 1:10
        f = x->evaluate_Bspline(N-1, x, T)
        # Integral should be 1
        if !(T <: BigFloat)
            I,e = QuadGK.quadgk(f, 0, N, rtol = tol)
            @test I≈T(1)
        end
        # Infinite summation of shifted versions is 1
        xx = linspace(T(N-1), T(N), S)[1:end-1]
        # xx = range(T(N-1), step=T(1)/T(S), length=S)
        g = zeros(T,length(xx))
        for k in 0:N-1
            g += map(x->f(x-k), xx)
        end
        @test g ≈ ones(T,length(g))  # (norm(g-1) < tol)
        # Two scale relation
        x = linspace(T(-1), T(N+1), S)
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
        for x in linspace(T(0),period,10)[1:end-1]
            @test (evaluate_periodic_Bspline(N, x, period, T) ≈ evaluate_Bspline(N, x, T))
            for k in -2:2
                @test (evaluate_periodic_Bspline(N, x, period, T) ≈ evaluate_periodic_Bspline(N, x+k*period, period, T))
            end
        end
        period = T(N+1)/T(3)
        for x in linspace(0,period,10)[1:end-1]
            for k in -2:2
                @test (evaluate_periodic_Bspline(N, x, period, T) ≈ evaluate_periodic_Bspline(N, x+k*period, period, T))
            end
        end
    end
end
function symmetricbsplinestest(T)
    K = 20
    for N in 0:4
        xs = linspace(T(0)+eps(T), T(1), K)
        for x in xs
            @test evaluate_symmetric_periodic_Bspline(N, x, T(10(N+1)), T) ≈ evaluate_symmetric_periodic_Bspline(N, -x, T(10(N+1)), T)
        end
    end
end

function test_spline_integration()
    T = Float64
    for N in 0:8
        @test squared_spline_integral(N) ≈ shifted_spline_integral(N,0)
        for t in 0:4
            f = x->evaluate_Bspline(N,x,T)
            @test abs(QuadGK.quadgk(x->f(x)*f(x+t),(-N-2:N+2)..., rtol=sqrt(eps(T)))[1]-shifted_spline_integral(N,t))<10*sqrt(eps(T))
        end
    end
end


function derivative_test(T)
    t = linspace(-10,10,100); t2 = (t[1:end-1]+t[2:end])/2
    for i in 2:10
        y1 = evaluate_Bspline.(i, t, T)
        y2 = diff_evaluate_Bspline.(i, t2, T)
        @test norm(diff(y1)./diff(t).-y2) < .1
    end
    for i in 2:10
        y1 = evaluate_periodic_Bspline.(i, t, 5, T)
        y2 = diff_evaluate_periodic_Bspline.(i, t2, 5, T)
        @test norm(diff(y1)./diff(t).-y2) < .1
    end

    for t in [linspace(0,1,100), linspace(-1,0,100), linspace(1,2,100)]
        t2 = (t[1:end-1]+t[2:end])/2
        y1 = evaluate_Bspline.(1, t, T)
        y2 = diff_evaluate_Bspline.(1, t2, T)
        @test norm(diff(y1)./diff(t).-y2) < .1
    end

    t0 = linspace(-1,7,100)
    for D in 0:6
        @test norm(diff_evaluate_Bspline.(Val{D}, Val{D+1}, t0, Float64)) == 0
    end

    t0 = -2.25:1/(1<<4):10
    t1 = (t0[1:end-1]+t0[2:end])/2
    t2 = (t1[1:end-1]+t1[2:end])/2
    t3 = (t2[1:end-1]+t2[2:end])/2
    t4 = (t3[1:end-1]+t3[2:end])/2
    t5 = (t4[1:end-1]+t4[2:end])/2

    for D in 3:6
        y0 = diff_evaluate_Bspline.(Val{D}, Val{0}, t0, Float64)
        y2 = diff_evaluate_Bspline.(Val{D}, Val{2}, t2, Float64)
        @test .5 > norm(diff(diff(y0)./diff(t0))./diff(t1) - y2)
    end

    for D in 4:6
        y0 = diff_evaluate_Bspline.(Val{D}, Val{0}, t0, Float64)
        y2 = diff_evaluate_Bspline.(Val{D}, Val{3}, t3, Float64)
        @test .5 > norm(  diff(diff(diff(y0)./diff(t0))./diff(t1))./diff(t2) - y2)
    end

    for D in 5:6
        y0 = diff_evaluate_Bspline.(Val{D}, Val{0}, t0, Float64)
        y2 = diff_evaluate_Bspline.(Val{D}, Val{4}, t4, Float64)
        @test .5 > norm(  diff(diff(diff(diff(y0)./diff(t0))./diff(t1))./diff(t2))./diff(t3) - y2)
    end

    for D in 6:8
        y0 = diff_evaluate_Bspline.(Val{D}, Val{0}, t0, Float64)
        y2 = diff_evaluate_Bspline.(Val{D}, Val{5}, t5, Float64)
        @test .5 > norm(  diff(diff(diff(diff(diff(y0)./diff(t0))./diff(t1))./diff(t2))./diff(t3))./diff(t4) - y2)
    end

    for D in 1:6
        @test Int.(diff_evaluate_Bspline.(Val{D}, Val{D}, -.5:D+2, T)) == vcat(0,[(-1)^k*binomial(D,k) for k in 0:D],0)
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
    @testset "$(rpad("symmetric B splines",P))"  begin
      symmetricbsplinestest(T)
    end
    @testset "$(rpad("derivative of B splines",P))"  begin
      derivative_test(T)
    end
end
@testset "$(rpad("integration of B splines",P))"  begin
  test_spline_integration()
end

if CREATE_README
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
