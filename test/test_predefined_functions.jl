@testset "Predefined Functions" begin
    @testset "lognormal" begin
        @test lognormal(1.0, 1.0, 1.0) == 0.0
        @test lognormal(2.0, 1.0, 1.0) == 0.5
        @test lognormal(1.0, 1.0, Inf) == 0
        @test lognormal(1.0, 1.0, 0.0) == Inf
    end

    @testset "logpoisson" begin
        @test logpoisson(10, 10.0) ≈ 0.0
        @test logpoisson(10, -1.0) == Inf
        @test logpoisson(10, 5.0, 5.0) ≈ 0.0
        @test logpoisson(10, 20.0) > 0
    end

    @testset "batlog" begin
        @test batlog(1.0) ≈ 0.0
        @test batlog(exp(1.0)) ≈ 1.0
        @test batlog(0.0) == -1e8
        @test batlog(-1.0) == -1e8
    end
end
