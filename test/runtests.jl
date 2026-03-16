using Batman
using Test

@testset "Batman.jl" begin
    include("poisson.jl")
    include("test_predefined_functions.jl")
    include("test_variable.jl")
    include("test_custom_model.jl")
end
