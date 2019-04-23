using CardinalBSplines, LinearAlgebra, Test, Sequences

@testset "Filter short cuts" begin
    for b in CardinalBSplines.IMPLEMENTED_FILTERS
        @test abs(b[0]) > abs(b[1])
    end
    for n in 0:3, m in 1:4, N in 10:3:20
        @test downsample(bnm(n,m,N)*b̃nm(n,m,N), m)[0:N-1]≈δ(0)[0:N-1]
        if m != 1
            @test downsample(bnm(n,m,N)*compact_dualperiodicbsplinesignal(n,m,N), m)[0:N-1]≈δ(0)[0:N-1]
        end
        @test abs(snm(n,m,N)[0]) > abs(snm(n,m,N)[1])
    end
end
