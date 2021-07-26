"""
Julia package designed to build statistical models with associated systematic
uncertainties and provide inference into particular parameters of interest.
Both Bayesian and Frequentist analysis are supported.
"""
module Batman 

using NLopt
using DataFrames
using DataStructures
using SpecialFunctions
using PyPlot
using Printf

# Logging
import Logging

function bat_format(level, _module, group, id, file, line)
  color, prefix, suffix = Logging.default_metafmt(level, _module, group, id, file, line)
  suffix = level == Logging.Debug ? "" : suffix
  return color, prefix, suffix
end

function bat_log(;kwargs...)
  kwargs = Dict(kwargs)
  debugLevel = haskey(ENV, "JULIA_DEBUG") ? Logging.Debug : Logging.Info
  debugLevel = get(kwargs, :level, debugLevel)
  Logging.ConsoleLogger( stdout, debugLevel; meta_formatter=bat_format )
end

Logging.global_logger( bat_log() )

# Primary definitions
include("predefined_functions.jl")
include("variable.jl")
include("likelihood.jl")
include("model.jl")
include("profile.jl")
include("plotter.jl")
include("show.jl")
# Kit definitions
include("modelkits/base.jl")
# Submodules
include("datastructures/DataStructures.jl")
include("toolbelt/Toolbelt.jl")
# Finally, macros
include("macros.jl")

# Batman exports based on a blacklist; where functions
# which begin with "_" are not exported.
const _EXCLUDE_SYMBOLS=[Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__, all=true)
  sym_string = string(sym)
  if sym in _EXCLUDE_SYMBOLS || startswith(sym_string, "_")
    continue
  end
  if !(Base.isidentifier(sym) || (startswith(sym_string, "@") &&
                                  Base.isidentifier(sym_string[2:end])))
    continue
  end
  @eval export $sym
end

end # End module Batman
