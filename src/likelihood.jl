import StatsBase
import Interpolations; itp=Interpolations
import KernelDensity

"""
    NLogPDF(f::String, p::Variable...)

Negative log PDF used when contructing the complete likelihood
function. The NLogPDF should have a unique name `f::String`.
"""
mutable struct NLogPDF 
  func::String
  params::Array{Variable}
  function NLogPDF(f::String, p::Variable...)
    new(f, [p...])
  end
end

"""
    NLogLikelihood(pf::Array{NLogPDF})

Given an array of `NLogPDF`, a complete negative log-likelihood
function is contructed which can then be passed to a minimizer.
"""
mutable struct NLogLikelihood
  objective::Function
  numparams::Int64
  variableList
  parameters
  lower_bounds
  upper_bounds
  objectiveString
  function NLogLikelihood(pf::Array{NLogPDF})
    var_list = []
    for p in pf
      for v in p.params
        push!(var_list, v)
      end
    end
    vv = unique(var_list)
    params = [v for v in vv if v.constant == false]
    lower_bounds = [v.bounds[1] for v in vv if v.constant == false]
    upper_bounds = [v.bounds[2] for v in vv if v.constant == false]
    nparams = length(params)
    npdf = length(pf)
    matches = []
    for (idx, p) in enumerate(pf)
      match = []
      for v in p.params
        loc = findlast(x->x==v, params)
        push!(match, loc)
      end
      push!(matches, match)
    end
    seval = ""
    for (idx, p) in enumerate(pf)
      a = matches[idx]
      name = pf[idx].func
      seval *= "$name("
      for (i,v) in zip(a, p.params)
        if v.constant
          val_use = v.init
          seval *= "$val_use,"
        else
          seval *= "x[$i],"
        end
      end
      seval *= ") +"
    end
    seval = seval[1:end-1]
    @debug "Likelihood function" seval
    ff = eval(Meta.parse("x->"*seval))
    new(ff, nparams, vv, params, lower_bounds, upper_bounds, seval)
  end
  function NLogLikelihood()
    new(x->x, 0, [], [], "")
  end
end

"""
    HistogramPDF(df::DataFrame, axis::Symbol; kwargs...)

Given a set of data `df`, along an axis `axis`, build a histogram
to be used as a PDF in a fit.

# Arguments
- `df::DataFrame`: Data to build the PDF
- `axis::Symbol`: Axis along which to slice the data
- `ndim::Integer=1`: Number of dimensions
- `bins=100`: Number of bins, or collection of bins.
- `interpolate=Interpolations.Constant()`: Constant or Linear interpolation
- `extrapolate=0`: Flat(), Linear(), or 0 extrapolation.
"""
function HistogramPDF(df::DataFrame, axis::Symbol; kwargs...)
  kwargs = Dict(kwargs)
  ndim = get(kwargs, :dimensions, 1)
  bins = get(kwargs, :bins, 100)
  itp_model = get(kwargs, :interpolate, itp.Constant())
  etp_model = get(kwargs, :extrapolate, 0)
  weight_symbol = get(kwargs, :weight, nothing)
  mass = weight_symbol != nothing ? StatsBase.weights(df[!, weight_symbol]) : StatsBase.weights(ones(size(df[!, axis])))
  h = StatsBase.fit( StatsBase.Histogram, df[!, axis], mass, bins; closed=:right )
  hx, hy = (h.edges[1][2:end] + h.edges[1][1:end-1])/2.0, h.weights
  #hx, hy = h.edges[1][2:end], h.weights

  # This works, but turning off for now
  return function_from_hist(h)

  ## SKIP BELOW
  hy = hy ./ sum(hy) ./ (hx[2]-hx[1])
  #hy = hy ./ sum(hy)
  #
  f = itp.extrapolate( 
        itp.interpolate((hx,), hy, itp.Gridded(itp_model)), etp_model
      )
  return x->f(x)
  # Now we need a more precise integral
  lowx = minimum(bins)
  hix = maximum(bins)
  bw = (hix - lowx) / length(bins) / 1_000_000.0
  integral_range = collect(lowx:bw:hix)
  integral = sum( f.(integral_range) )*bw
  x->f(x)/integral
end

function function_from_hist(h::StatsBase.Histogram)
  # Normalize
  hx, hy = h.edges, h.weights
  hx = hx[1]
  hy = hy ./ sum( hy ) ./ (hx[2] - hx[1])
  push!(hy, 0)
  newfunc = function (x)
    lasty = 0
    for (testx, testy) in zip(hx, hy)
      ## JUST CHANGED TO >=?
      if testx > x
        return lasty
      end
      lasty = testy
    end
    0
  end
  return x -> newfunc.(x)
  ## We should test the function real quick I suppose
  #bw = (hx[2] - hx[1]) / 1_000.0
  #integral_range = collect(hx[1]:bw:hx[end])
  #integral = sum( newfunc.(integral_range) )*bw
  #@show "New Integral:: " integral
  #newnewfunc = x->newfunc.(x) / integral
  #bw = (hx[2] - hx[1]) / 10_000.0
  #integral_range = collect(hx[1]:bw:hx[end])
  #integral = sum( newnewfunc.(integral_range) )*bw
  #@show "New New Integral:: " integral
  #return x -> newnewfunc.(x)
end

## Okay, lets throw this in a temp funciton
#function HistogramPDF(df::DataFrame, axis::Symbol; kwargs...)
#  a = rand()*10.0
#  b = rand()*3.0
#  #(x)->(a.*x .+ b)./(a+b)
#  norm = (cos(b)-cos(a+b)+a)/(2*a)
#  (x)-> (sin.(a.*x .+ b) .+ 1) ./ 2 ./ norm
#end

"""
    KdePDF(df::DataFrame, axis::Symbol; kwargs...)

Given a set of data `df`, along an axis `axis`, produce a kernel
density estimate to use as a PDF in a fit.

# Arguments
- `df::DataFrame`: Data to build the PDF
- `axis::Symbol`: Axis along which to slice the data
- `ndim::Integer=1`: Number of dimensions
"""
function KdePDF(df::DataFrame, axis::Symbol; kwargs...)
  kwargs = Dict(kwargs)
  # :dimensions not used currently
  ndim = get(kwargs, :dimensions, 1)

  k = KernelDensity.kde(df[!, axis])
  imd = KernelDensity.InterpKDE(k)
  p(x) = KernelDensity.pdf(imd, x)
  return p
end
