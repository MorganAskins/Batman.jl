import Distributions: Poisson

struct SingleSpectra
    f
    p::Parameter
    o::Observable
    df::Symbol
end

"""
    Spectralfit()

True multidimensional fit in N-Dimensions, derived from Monte Carlo histograms
of data.

# Example
```julia
m = Spectralfit()
```
"""
mutable struct Spectralfit <: Model
    params::Vector
    pdflist::Vector{NLogPDF}
    nll::NLogLikelihood
    spectra::Vector{SingleSpectra}
    datasets::Vector{Symbol}
    observables::Vector{Observable}
    function Spectralfit(observables::Symbol...)
        new([], NLogPDF[], NLogLikelihood(),
            SingleSpectra[], Symbol[], Observable[])
    end
end

struct SpectralPDF
    p::Parameter
    f::String
end

"""
    generate_mock_dataset(m::Spectralfit; time=1.0)

Generate a mock dataset for the model using Poisson sampling.
"""
function generate_mock_dataset(m::Spectralfit; kwargs...)
    kwargs = Dict(kwargs)
    time = get(kwargs, :time, 1.0)
    for d in m.datasets
        df = eval(d)
        delete!(df, 1:size(df, 1))
    end
    for p in m.params
        p.fit = rand(Poisson(p.init))
    end
    for sp in m.spectra
        append!(sp.df[!, sp.o.sname],
                randomFromSpectrum(sp.f; count=sp.p.fit, min=sp.o.min, max=sp.o.max))
    end
end

mutable struct EventCounter
    name::String
    count::Int
    spectra::Vector{SingleSpectra}
    obs::Dict{Symbol,Vector{Float64}}
    function EventCounter(name::String, count::Int)
        new(name, count, SingleSpectra[], Dict{Symbol,Vector{Float64}}())
    end
end

function build_dataframe(ev::EventCounter)
    df = DataFrame()
    for (k, v) in ev.obs
        df[!, k] = v
    end
    df[!, :meta] = repeat([ev.name], ev.count)
    return df
end

function compare_and_add!(df1::DataFrame, df2)
    names1 = names(df1)
    names2 = names(df2)
    for n1 in names1
        if !(n1 in names2)
            df2[!, n1] = typeof(df1[!, n1])(undef, size(df2, 1))
        end
    end
    for n2 in names2
        if !(n2 in names1)
            df1[!, n2] = typeof(df2[!, n2])(undef, size(df1, 1))
        end
    end
    append!(df1, df2)
end

function genmo(m::Spectralfit; kwargs...)
    kwargs = Dict(kwargs)
    time = get(kwargs, :time, 1.0)
    verbose = get(kwargs, :verbose, false)
    countdict = get(kwargs, :counts, Dict())

    for d in m.datasets
        df = eval(d)
        delete!(df, 1:size(df, 1))
        @debug "Deleting from" d
    end
    dictDPCount = Dict{Symbol,Dict{String,EventCounter}}(sp.df => Dict() for sp in m.spectra)
    for sp in m.spectra
        param = sp.p
        name = param.name
        count_init = get(countdict, name, param.init)
        count::Int = rand(Poisson(param.efficiencies[sp.df] * time * count_init))
        if verbose
            println("Name: ", param.name, " Eff: ", param.efficiencies,
                    " Init: ", param.init, " time ", time, " count ", count)
        end
        if !haskey(dictDPCount[sp.df], name)
            dictDPCount[sp.df][name] = EventCounter(name, count)
        end
        push!(dictDPCount[sp.df][name].spectra, sp)
    end

    for (dset, evdict) in dictDPCount
        for (_, ev) in evdict
            count = ev.count
            for sp in ev.spectra
                newdata = randomFromSpectrum(sp.f; count=count, min=sp.o.min, max=sp.o.max)
                ev.obs[sp.o.sname] = newdata
            end
            dfnew = build_dataframe(ev)
            compare_and_add!(getfield(Batman, dset), dfnew)
        end
    end
    updateDataFunctions!(m)
    return dictDPCount
end

function updateDataFunctions!(m::Spectralfit)
    for sp in m.spectra
        pname = sp.p.name
        obs = sp.o.name
        fSymbol = Symbol(pname * "_" * obs * "_" * String(sp.df))
        add_array(fSymbol, sp.f(eval(sp.df)[!, sp.o.sname]))
    end
end

function constructPDF!(m::Spectralfit, p::Parameter, shapes, observables, ds::Symbol)
    if !(ds in m.datasets)
        push!(m.datasets, ds)
    end
    seval = ""
    for (shp, obs) in zip(shapes, observables)
        fSymbol = Symbol(p.name * "_" * String(obs) * "_" * String(ds))
        add_array(fSymbol, shp(eval(ds)[!, obs]))
        push!(m.spectra, SingleSpectra(shp, p, observableFromSymbol(m, obs), ds))
        seval *= "$fSymbol.*"
    end
    seval = seval[1:end-2]
    @debug seval
    unique_name = "spectralPDF_" * p.name * "_" * String(ds)
    add_function(unique_name, eval(Meta.parse("()->" * seval)))
    SpectralPDF(p, unique_name)
end

function combinePDFs!(m::Spectralfit, spectralpdfs::Vector{SpectralPDF},
                      ds::Symbol; kwargs...)
    kwargs = Dict(kwargs)
    livetime = get(kwargs, :livetime, nothing)
    eff = get(kwargs, :efficiency, Vector{Nothing}(undef, length(spectralpdfs)))

    seval = "-sum(batlog.("
    fcall = "("
    extended = "+sum("
    plist = []
    if livetime !== nothing
        push!(plist, livetime)
        tname = livetime.name
        fcall *= "$tname,"
    end
    for (spdf, e) in zip(spectralpdfs, eff)
        p = spdf.p
        f = spdf.f
        vname = p.name
        if e !== nothing
            ename = e.name
            seval *= "$ename*"
            fcall *= "$ename,"
            extended *= "$ename*"
            push!(plist, e)
        end
        push!(plist, p)
        if livetime !== nothing
            tname = livetime.name
            seval *= "$tname*"
            extended *= "$tname*"
        end
        vinit = p.init
        seval *= "($vname)*$vinit*$f().+"
        fcall *= "$vname,"
        extended *= "($vname)*$vinit+"
    end
    seval = seval[1:end-2]
    seval *= "))"
    fcall = fcall[1:end-1]
    fcall *= ")->"
    extended = extended[1:end-1]
    extended *= ")" * "+size($ds,1)*log(size($ds,1))-size($ds,1)"
    unique_name = "fullSpectralPDF_" * String(ds)
    @debug fcall * seval * extended
    add_function(unique_name, eval(Meta.parse(fcall * seval * extended)))
    push!(m.pdflist, NLogPDF(unique_name, plist...))
end

function add_parameter!(m::Spectralfit, name::String, init; σ=Inf, kwargs...)
    kwargs = Dict(kwargs)
    dsEfficiency = get(kwargs, :dsEfficiency, (Symbol(name) => 1.0))
    bounds = get(kwargs, :bounds, (-Inf, Inf))
    p = getparam(m, name)
    if p !== nothing
        @warn name * " already exists in model."
        push!(p.efficiencies, dsEfficiency)
        return p
    end
    @debug "Init" init
    p = Parameter(name; init=init, bounds=bounds)
    push!(p.efficiencies, dsEfficiency)
    p.samplewidth = σ
    push!(m.params, p)
    if σ != Inf
        a = Constant(name * "_mu", 1.0)
        b = Constant(name * "_sig", σ)
        nlp = NLogPDF("lognormal", p, a, b)
        push!(m.pdflist, nlp)
    end
    return p
end

function add_constant!(m::Spectralfit, name::String, init)
    k = getparam(m, name)
    if k !== nothing
        @warn name * " already exists in model."
        return k
    end
    Constant(name, init)
end

function add_observable!(m::Spectralfit, name::Symbol, min::Float64, max::Float64)
    push!(m.observables, Observable(String(name), name, min, max))
end

function observableFromSymbol(m::Spectralfit, name::Symbol)
    for obs in m.observables
        if obs.sname == name
            return obs
        end
    end
    return nothing
end

function minimize!(m::Spectralfit; options=Dict())
    likelihood = NLogLikelihood(m.pdflist)
    add_likelihood!(m, likelihood)
    optimize_model!(m, likelihood; options=options)
end

"""
    randomFromSpectrum(f; min=0.0, max=1.0, top=10.0, count=1)

Return a `Vector{Float64}` of random elements from a 1D PDF given by `f`
between the values of `min` and `max` with a maximum probability of `top`.
"""
function randomFromSpectrum(f; min=0.0, max=1.0, top=10.0, count=1)
    xArr = Float64[]
    for _ in 1:count
        while true
            x = rand() * (max - min) + min
            y = rand() * top
            if y < f(x)
                push!(xArr, x)
                break
            end
        end
    end
    return xArr
end
