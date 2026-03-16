"""
    rootreader(fileName::String, treeName::String; kwargs...)

Read a single ROOT TTree from a single TFile and return
the data as a DataFrame; preserving the names of the variables
and order of the events.

# Keywords
- `cuts::Vector{Expr}`: Slice rows from the dataframe given cuts on column values.
- `observables::Vector{Symbol}`: Select to keep specific columns.
- `meta::Dict`: Internal information stored in this dict if passed in.

# Example
```julia
meta = Dict()
rootreader("/path/to/file.root", "treename";
    cuts=[:(energy .> 5.5), :(position .< 1.5)],
    observables=[:energy, :position],
    meta=meta
)
```

!!! note
    Requires UpROOT.jl to be installed separately.
"""
function rootreader(fileName::String, treeName::String; kwargs...)
    kwargs = Dict(kwargs)
    cols = get(kwargs, :cuts, Expr[])
    rows = get(kwargs, :observables, nothing)
    meta = get(kwargs, :meta, Dict())

    uproot = Base.require(Main, :UpROOT).uproot
    tfile    = uproot.open(fileName)
    ttree    = tfile.get(treeName)
    branches = rows === nothing ? ttree.keys() : [String(r) for r in rows]
    events   = ttree.arrays(branches)
    df = DataFrame(events)
    meta["entries"] = size(df, 1)
    for c in cols
        df = rowselect(df, c)
    end
    return df
end
