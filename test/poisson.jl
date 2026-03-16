@testset "CountingExperiment" begin
    m = CountingExperiment()
    add_component!(m, "Signal", 20.0)
    add_component!(m, "Bkg 1", 30.0; σ=0.1)
    add_component!(m, "Bkg 2", 30.0; σ=45.0)
    add_component!(m, "Bkg 3", 12.0; σ=12.0)

    set_counts!(m, 100)

    results = minimize!(m)

    compute_profiled_uncertainties!(results; σ=1)

    signal = getparam(m, "Signal")

    @test signal.fit ≈ 18 rtol = 0.5
    @test signal.low < signal.fit
    @test signal.high > signal.fit
    @test results.min_objective isa Float64
    @test results.iterations > 0
end
