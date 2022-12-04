"""
    CountingExperiment()

Building a Poisson counting model of n>1 components, each
added individual, including uncertainties.

# Example
```jldoctest
julia> m = CountingExperiment()

julia> add_component!(m, "Signal", 20.0)
julia> add_component!(m, "Background", 30.0; σ=5.0)

julia> set_counts!(m, 55)
julia> results = minimize!(m)
```
"""
mutable struct CountingExperiment <: Model
  params
  pdflist
  nll::NLogLikelihood
  counts::Constant
  function CountingExperiment()
    new( [], [], NLogLikelihood(),
        Constant("counts", 0.0) )
  end
end

function add_component!( m::CountingExperiment, name::String, μ; σ=Inf )
  p = Parameter( name; init=μ, bounds=(0.0, Inf) )
  push!(m.params, p)
  if σ != Inf
    a = Constant( name*"_mu", μ )
    b = Constant( name*"_sig", σ )
    nlp = NLogPDF( "lognormal", p, a, b )
    push!(m.pdflist, nlp)
  end
end

function set_counts!( m::CountingExperiment, counts::Number )
  m.counts.init = Float64(counts)
end

function minimize!( m::CountingExperiment )
  poisson_pdf = NLogPDF("logpoisson", m.counts, (m.params)...)
  likelihood = NLogLikelihood([m.pdflist..., poisson_pdf])
  # Fitter options
  add_likelihood!( m, likelihood )
  optimize_model!( m, likelihood )
  #optimize_model!( m, likelihood; 
  #                 lower_bounds=[0.0 for a in 1:length(m.params)] )
end
