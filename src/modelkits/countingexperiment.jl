"""
    CountingExperiment()

Build a Poisson counting model of n>1 components, each
added individually, including uncertainties.

# Example
```julia
m = CountingExperiment()
add_component!(m, "Signal", 20.0)
add_component!(m, "Background", 30.0; σ=5.0)
set_counts!(m, 55)
results = minimize!(m)
```
"""
mutable struct CountingExperiment <: Model
    params::Vector
    pdflist::Vector{NLogPDF}
    nll::NLogLikelihood
    counts::Constant
    function CountingExperiment()
        new([], NLogPDF[], NLogLikelihood(), Constant("counts", 0.0))
    end
end

function add_component!(m::CountingExperiment, name::String, μ; σ=Inf)
    p = Parameter(name; init=μ, bounds=(0.0, Inf))
    push!(m.params, p)
    if σ != Inf
        a = Constant(name * "_mu", μ)
        b = Constant(name * "_sig", σ)
        nlp = NLogPDF("lognormal", p, a, b)
        push!(m.pdflist, nlp)
    end
end

function set_counts!(m::CountingExperiment, counts::Number)
    m.counts.init = Float64(counts)
end

function minimize!(m::CountingExperiment)
    poisson_pdf = NLogPDF("logpoisson", m.counts, (m.params)...)
    likelihood = NLogLikelihood([m.pdflist..., poisson_pdf])
    add_likelihood!(m, likelihood)
    optimize_model!(m, likelihood)
end
