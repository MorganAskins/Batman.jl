@testset "CustomModel" begin
    m = CustomModel()
    @test isempty(m.params)
    @test isempty(m.pdflist)

    p1 = Parameter("x"; init=5.0, bounds=(0.0, 100.0))
    mu = Constant("mu", 5.0)
    sig = Constant("sig", 1.0)

    add_parameters!(m, [p1])
    @test length(m.params) == 1

    nlp = NLogPDF("lognormal", p1, mu, sig)
    add_nlogpdfs!(m, [nlp])
    @test length(m.pdflist) == 1

    results = minimize!(m)
    @test results.min_parameters[1] ≈ 5.0 rtol = 0.1
end
