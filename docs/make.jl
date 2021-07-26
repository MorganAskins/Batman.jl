push!(LOAD_PATH, "../src/")
using Documenter, Batman 

makedocs(
  sitename="Batman.jl",
  modules=[Batman],
  pages = Any[
    "Home" => "index.md",
    "Library" => map(s -> "lib/$(s)", sort(readdir(joinpath(@__DIR__, "src/lib"))))
  ],
  format = Documenter.HTML( prettyurls = get(ENV, "CI", nothing) == "true" )
)

deploydocs(
  repo = "github.com/MorganAskins/Batman.jl.git",
)
