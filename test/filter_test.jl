using CardinalBSplines, LinearAlgebra, Test

@testset "Filter short cuts" begin
    for b in CardinalBSplines.IMPLEMENTED_FILTERS
        @test abs(b[0]) > abs(b[1])
    end
end
