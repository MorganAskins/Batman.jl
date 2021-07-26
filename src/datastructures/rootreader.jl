import UpROOT: UpROOT.uproot

"""
rootreader(fileName::String, treeName::String)

Read a single ROOT TTree from a single TFile and return
the data as a DataFrame; preserving the names of the variables
and order of the events.

# Arguments
- `cuts::Array{Expr}`: Slice rows from the dataframe given cuts on column values.
- `observables::Array{Symbol}`: Select to keep specific columns.
- `meta::Dict`: Internal information stored in this dict if passed in.

# Example
```julia
julia> meta = Dict()
julia> rootreader("/path/to/file.root", "treename";
         cuts=[:(energy .> 5.5), :(position .< 1.5)],
         observables=[:energy, :position],
         meta=meta
       );
```
"""
function rootreader(fileName::String, treeName::String; kwargs...)
  kwargs = Dict(kwargs)
  # Options: Apply cuts; eg
  # cuts = ( (:energy >. 5.0) & (:posr3 <. 1.0) )
  cols = get(kwargs, :cuts, Array{Expr}(undef, 0) )
  rows = get(kwargs, :observables, nothing)
  meta = get(kwargs, :meta, Dict())

  tfile    = uproot.open(fileName)
  ttree    = tfile.get(treeName)
  branches = rows == nothing ? ttree.keys() : [String(r) for r in rows]
  events   = ttree.arrays(branches)
  df = DataFrame(events)
  meta["entries"] = size(df, 1)
  for c in cols
    df = rowselect(df, c)
  end
  df
end
