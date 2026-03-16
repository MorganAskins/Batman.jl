using Documenter
using Batman

makedocs(;
    sitename="Batman.jl",
    modules=[Batman],
    pages=[
        "Home" => "index.md",
        "Library" => map(s -> "lib/$(s)", sort(readdir(joinpath(@__DIR__, "src/lib")))),
    ],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://MorganAskins.github.io/Batman.jl",
    ),
    warnonly=true,
)

deploydocs(;
    repo="github.com/MorganAskins/Batman.jl.git",
    devbranch="master",
    push_preview=true,
)
