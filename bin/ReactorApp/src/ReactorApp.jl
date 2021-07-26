#!/bin/env julia

module ReactorApp

push!(LOAD_PATH, "../src")
using Batman
using JSON
using ArgParse
using PyPlot
using SpecialFunctions

mutable struct reactor_storage
  number_components::Int64
  parameter_names::Array{String}
  parameter_inits::Array{Float64}
  parameter_errors::Array{Float64}
  options::Dict{String,Any}
end

function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table! s begin
    "--input", "-i"
      help="Input configuration file"
      required=true
    "--verbose", "-v"
      help="Increase verbosity"
      action=:store_true
  end
  return parse_args(s)
end

function parse_json(args)
  ifile = args["input"]
  data = JSON.parsefile(ifile)
  config_options = data["options"]
  par_names::Array{String} = Array{String}(undef, 0)
  par_inits::Array{Float64} = Array{Float64}(undef, 0)
  par_errors::Array{Float64} = Array{Float64}(undef, 0)
  for (k,v) in data["components"]
    push!(par_names, k)
    push!(par_inits, v["rate"])
    push!(par_errors, v["error"])
  end
  reactor_storage( length(data["components"]),
                   par_names,
                   par_inits,
                   par_errors,
                   config_options )
end

function sensitivity(name, results;
                     mode="FC", prior=nothing)
  comp = getparam(results.model, name)
  x, y = comp.likelihood_x, comp.likelihood_y
  if prior != nothing
    y = y .* prior.(x)
  end
  xselect = (x .>= 0)
  x = x[xselect]
  y = y[xselect]
  y = y / sum(y)
  if mode == "FC"
    cy = sortperm(y; rev=true)
    cum_y = cumsum(y[cy])
    cum_x = x[cy]
    first_min = cum_y[ cum_x .== x[1] ]
  end
  return minimum([first_min[1], 1])
end

"""
This script is intended to run as an executable which reads
in a json configuration file and outputs a series of objects
such as plots, tables, and a summary report.
Todo:
> 
"""
function julia_main()
  parsed_args    = parse_commandline()
  reactor_config = parse_json(parsed_args)
  options = reactor_config.options
  ## Gather options
  nsigma = get(options, "nsigma", 6)
  poi = get(options, "parameters_of_interest", "signal")
  days = get(options, "days", 100)
  hs = x -> x >= 0 ? 1 : 0
  results_set = []

  days = range(1, length=days, stop=days)
  for d in days
    ## Build a Poisson Model
    m = CountingExperiment()
    n_components = reactor_config.number_components
    asimov = 0
    for i in collect(1:n_components)
      error = reactor_config.parameter_errors[i]
      error = error > 0 ? error : Inf
      add_component!(m, reactor_config.parameter_names[i], 
                        reactor_config.parameter_inits[i]*d;
                        σ=error*d)
      asimov += reactor_config.parameter_inits[i]
    end
    set_counts!(m, asimov*d)

    results = minimize!(m)
    profile!(poi, results; prior=hs, step=0.1, stop=20)
    uncertainty!(poi, results, σ=1)
    push!(results_set, results)
  end

  sens = []
  for rs in results_set
    push!(sens, erfinv(sensitivity(poi, rs))*2^0.5)
  end
  @show length(days)
  @show length(sens)
  plt.plot(days, sens)
  plt.savefig("results.svg")
  return 0
end

end # module
