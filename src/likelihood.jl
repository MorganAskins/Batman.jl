import StatsBase
import Interpolations; const itp = Interpolations
import KernelDensity

"""
    NLogPDF(f::String, p::Variable...)

Negative log PDF used when constructing the complete likelihood
function. The NLogPDF should have a unique name `f::String`.
"""
mutable struct NLogPDF
    func::String
    params::Vector{Variable}
    function NLogPDF(f::String, p::Variable...)
        new(f, [p...])
    end
end

"""
    NLogLikelihood(pf::Vector{NLogPDF})

Given a vector of `NLogPDF`, a complete negative log-likelihood
function is constructed which can then be passed to a minimizer.
"""
mutable struct NLogLikelihood
    objective::Function
    numparams::Int
    variableList::Vector
    parameters::Vector
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    objectiveString::String
    function NLogLikelihood(pf::Vector{NLogPDF})
        var_list = Variable[]
        for p in pf
            for v in p.params
                push!(var_list, v)
            end
        end
        vv = unique(var_list)
        params = [v for v in vv if !v.constant]
        lower_bounds = [v.bounds[1] for v in vv if !v.constant]
        upper_bounds = [v.bounds[2] for v in vv if !v.constant]
        nparams = length(params)
        matches = Vector{Vector{Union{Int,Nothing}}}()
        for p in pf
            match = Union{Int,Nothing}[]
            for v in p.params
                loc = findlast(x -> x == v, params)
                push!(match, loc)
            end
            push!(matches, match)
        end
        seval = ""
        for (idx, p) in enumerate(pf)
            a = matches[idx]
            name = pf[idx].func
            seval *= "$name("
            for (i, v) in zip(a, p.params)
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
        ff = eval(Meta.parse("x->" * seval))
        new(ff, nparams, vv, params, lower_bounds, upper_bounds, seval)
    end
    function NLogLikelihood()
        new(x -> x, 0, [], [], Float64[], Float64[], "")
    end
end

function FunctionPDF(fun; kwargs...)
    kwargs = Dict(kwargs)
    return x -> fun.(x...)
end

"""
    MultiHistogramPDF(df::DataFrame, axes::Vector{Symbol}; kwargs...)

Given a set of data `df`, along axes `axes`, build a histogram
to be used as a PDF in a fit.

# Keywords
- `dimensions::Int=1`: Number of dimensions
- `bins=100`: Number of bins, or collection of bins.
- `interpolate=Interpolations.Constant()`: Constant or Linear interpolation
- `extrapolate=0`: Flat(), Linear(), or 0 extrapolation.
- `weight=nothing`: Column symbol for weights.
- `probability="Histogram"`: PDF estimation method ("Histogram" or "KDE").
"""
function MultiHistogramPDF(df::DataFrame, axes::Vector{Symbol}; kwargs...)
    kwargs = Dict(kwargs)
    bins = Tuple(get(kwargs, :bins, 100))
    itp_model = get(kwargs, :interpolate, itp.Constant())
    etp_model = get(kwargs, :extrapolate, 0)
    weight_symbol = get(kwargs, :weight, nothing)
    probability = get(kwargs, :probability, "Histogram")

    z = Tuple([df[!, s] for s in axes])
    mass = weight_symbol !== nothing ? StatsBase.weights(df[!, weight_symbol]) : StatsBase.weights(ones(size(df[!, axes[1]])))
    if probability == "KDE"
        k = KernelDensity.kde(z)
        imd = KernelDensity.InterpKDE(k)
        p(x) = KernelDensity.pdf(imd, x)
        return x -> p.(x...)
    else
        h = StatsBase.fit(StatsBase.Histogram, z, mass, bins)
        hxs = Tuple([(j[2:end] + j[1:end-1]) / 2.0 for j in h.edges])
        volume = prod([(j[2] - j[1]) for j in hxs])
        hy = h.weights
        hy = hy ./ sum(hy) ./ volume
        f = itp.extrapolate(itp.interpolate(hxs, hy, itp.Gridded(itp_model)), etp_model)
        return x -> f.(x...)
    end
end

"""
    HistogramPDF(df::DataFrame, axis::Symbol; kwargs...)

Given a set of data `df`, along an axis `axis`, build a histogram
to be used as a PDF in a fit.

# Keywords
- `dimensions::Int=1`: Number of dimensions
- `bins=100`: Number of bins, or collection of bins.
- `interpolate=Interpolations.Constant()`: Constant or Linear interpolation
- `extrapolate=0`: Flat(), Linear(), or 0 extrapolation.
- `weight=nothing`: Column symbol for weights.
"""
function HistogramPDF(df::DataFrame, axis::Symbol; kwargs...)
    kwargs = Dict(kwargs)
    bins = get(kwargs, :bins, 100)
    itp_model = get(kwargs, :interpolate, itp.Constant())
    etp_model = get(kwargs, :extrapolate, 0)
    weight_symbol = get(kwargs, :weight, nothing)
    mass = weight_symbol !== nothing ? StatsBase.weights(df[!, weight_symbol]) : StatsBase.weights(ones(size(df[!, axis])))
    h = StatsBase.fit(StatsBase.Histogram, df[!, axis], mass, bins)
    hx, hy = (h.edges[1][2:end] + h.edges[1][1:end-1]) / 2.0, h.weights
    hy = hy ./ sum(hy) ./ (hx[2] - hx[1])
    f = itp.extrapolate(
        itp.interpolate((hx,), hy, itp.Gridded(itp_model)), etp_model
    )
    return x -> f(x)
end

"""
    KdePDF(df::DataFrame, axis::Symbol; kwargs...)

Given a set of data `df`, along an axis `axis`, produce a kernel
density estimate to use as a PDF in a fit.

# Keywords
- `dimensions::Int=1`: Number of dimensions (currently unused).
"""
function KdePDF(df::DataFrame, axis::Symbol; kwargs...)
    k = KernelDensity.kde(df[!, axis])
    imd = KernelDensity.InterpKDE(k)
    p(x) = KernelDensity.pdf(imd, x)
    return p
end
