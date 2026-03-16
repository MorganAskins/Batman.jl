@testset "Variables" begin
    @testset "Parameter" begin
        p = Parameter("test"; init=5.0, bounds=(0.0, 10.0))
        @test p.name == "test"
        @test p.init == 5.0
        @test p.bounds == (0.0, 10.0)
        @test p.constant == false
        @test isempty(p.likelihood_x)
        @test isempty(p.likelihood_y)
    end

    @testset "Constant" begin
        c = Constant("fixed", 3.14)
        @test c.name == "fixed"
        @test c.init == 3.14
        @test c.constant == true
    end

    @testset "Dataset" begin
        d = Dataset("mydata")
        @test d.name == "mydata"
        @test d.constant == true
    end

    @testset "Observable" begin
        o = Observable("energy", :energy, 0.0, 10.0)
        @test o.name == "energy"
        @test o.sname == :energy
        @test o.min == 0.0
        @test o.max == 10.0
    end
end
