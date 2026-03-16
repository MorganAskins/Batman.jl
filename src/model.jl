"""
    Model

Abstract type for all Batman models. Models hold the parameters of interest
and the associated likelihood/cost functions used in optimization.
"""
abstract type Model end

"""
    Results

Created upon optimizing a model, contains the model and all of
its parameters as well as the optimized results found by
minimizing the cost function.
"""
mutable struct Results
    opt::Opt
    model::Model
    min_objective::Float64
    min_parameters::Vector{Float64}
    iterations::Int
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
end

"""
    optimize_model!(m::Model, nll::NLogLikelihood; options=Dict())

Minimize the objective function from the given `NLogLikelihood` using NLopt.
Returns a `Results` object.

# Options (via `options` Dict)
- `"optimizer"`: One of "sbplx", "cobyla", "bobyqa", "praxis", "nelder"
- `"p0"`: Initial parameter values
- `"trials"`: Number of random restarts (default 1)
- `"ftol_abs"`, `"ftol_rel"`, `"xtol_rel"`, `"xtol_abs"`: Tolerances
- `"maxeval"`, `"maxtime"`: Stopping criteria
- `"initial_step"`: Initial step sizes
- `"stopval"`: Stop when objective reaches this value
"""
function optimize_model!(m::Model, nll::NLogLikelihood;
                         lower_bounds=nothing, upper_bounds=nothing,
                         options=Dict())
    optimizer_selector = Dict(
        "sbplx"  => :LN_SBPLX,
        "cobyla" => :LN_COBYLA,
        "bobyqa" => :LN_BOBYQA,
        "praxis" => :LN_PRAXIS,
        "nelder" => :LN_NELDERMEAD
    )
    opt_string = get(options, "optimizer", "sbplx")
    @debug "Using optimizer" opt_string
    opt_use = optimizer_selector[opt_string]
    function objective(x::Vector, grad::Vector)
        if length(grad) > 0
            grad .= 2x
        end
        nll.objective(x)
    end
    opt = Opt(opt_use, nll.numparams)
    p0 = get(options, "p0", [v.init for v in nll.parameters])
    trials = get(options, "trials", 1)

    opt.ftol_abs     = get(options, "ftol_abs", 1e-6)
    opt.ftol_rel     = get(options, "ftol_rel", 0)
    opt.xtol_rel     = get(options, "xtol_rel", 0)
    opt.xtol_abs     = get(options, "xtol_abs", 0)
    opt.maxeval      = get(options, "maxeval", 0)
    opt.maxtime      = get(options, "maxtime", 0)
    initial_step     = get(options, "initial_step", p0 ./ 10.0)
    initial_step[initial_step .== 0] .= 0.1
    opt.initial_step = initial_step
    opt.stopval      = get(options, "stopval", -Inf)
    opt.lower_bounds = nll.lower_bounds
    opt.upper_bounds = nll.upper_bounds

    opt.min_objective = objective

    best_minf = 1e9
    if trials > 1
        minf_set = Float64[]
        p0_set = Vector{Float64}[]
        for _ in 1:trials
            p0 = Float64[]
            for v in nll.parameters
                if v.samplewidth != Inf
                    push!(p0, randn() * v.samplewidth + v.init)
                else
                    push!(p0, v.init)
                end
            end
            (minf, _, _) = optimize!(opt, p0)
            push!(p0_set, p0)
            push!(minf_set, minf)
        end
        p0 = p0_set[minf_set .== minimum(minf_set)][1]
        best_minf = minimum(minf_set)
    end

    initial_step[initial_step .== 0] .= 0.1
    opt.initial_step = initial_step
    (minf, minx, _) = optimize!(opt, p0)

    if minf > best_minf
        @debug "Trial minimum was better" best_minf minf
    end

    numevals = opt.numevals
    for (i, p) in enumerate(m.params)
        p.fit = minx[i]
    end
    Results(opt, m, copy(minf), copy(minx), numevals,
            copy(opt.lower_bounds), copy(opt.upper_bounds))
end

function getparam(m::Model, name::String)
    for p in m.params
        if p.name == name
            return p
        end
    end
    return nothing
end

function add_likelihood!(m::Model, nll::NLogLikelihood)
    m.nll = nll
    m.params = nll.parameters
end
