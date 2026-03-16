"""
    profile!(name, results; kwargs...)

Compute the profiled likelihood over a single parameter of
interest given by `name`.

# Keywords
- `stop::Float64=5.0`: Maximum ΔLog(Likelihood)
- `step::Float64=1.0`: Incremental step of the parameter of interest
- `prior=nothing`: Prior probability of the parameter of interest
- `verbose::Bool=false`: Print diagnostic information
- `maxsteps::Int=1000`: Maximum number of profile steps
"""
function profile!(name, results; kwargs...)
    kwargs = Dict(kwargs)
    stop      = get(kwargs, :stop, 5.0)
    step      = get(kwargs, :step, 1.0)
    prior     = get(kwargs, :prior, nothing)
    verbose   = get(kwargs, :verbose, false)
    max_steps = get(kwargs, :maxsteps, 1000)

    fitterOptions = get(kwargs, :fitteroptions, nothing)
    if fitterOptions !== nothing
        for (k, v) in fitterOptions
            setproperty!(results.opt, k, v)
        end
    end

    idx = 0
    for (i, p) in enumerate(results.model.params)
        if p.name == name
            idx = i
        end
    end
    if idx == 0
        @warn "Parameter not found" name
        return Float64[], Float64[]
    end

    verb_count = 0
    verb_steps = Float64[]
    param = results.model.params[idx]

    # Scan left
    x_left = Float64[]
    nll_left = Float64[]
    nlow = copy(results.lower_bounds)
    nhigh = copy(results.upper_bounds)
    p0 = copy(results.min_parameters)
    val = copy(results.min_objective)
    xeval = p0[idx]
    count = 0
    while (val - results.min_objective) < stop
        nlow[idx] = xeval
        nhigh[idx] = xeval
        results.opt.lower_bounds = nlow
        results.opt.upper_bounds = nhigh
        p0[idx] = xeval
        init_step = get(kwargs, :init_step, p0 ./ 100.0)
        init_step[idx] = xeval * 1e-6
        init_step[init_step .== 0] .= 0.1
        results.opt.initial_step = init_step
        val, _, _ = optimize!(results.opt, p0)
        likelihood_value = exp(-(val - results.min_objective))
        push!(verb_steps, results.opt.numevals)
        if (likelihood_value == Inf) || (likelihood_value < 0)
            break
        end
        push!(nll_left, likelihood_value)
        push!(x_left, xeval)
        xeval -= step
        count += 1
        if count > max_steps
            @warn "Profiling hit limit" max_steps
            break
        end
    end
    verb_count += count

    # Scan right
    x_right = Float64[]
    nll_right = Float64[]
    nlow = copy(results.lower_bounds)
    nhigh = copy(results.upper_bounds)
    p0 = copy(results.min_parameters)
    val = copy(results.min_objective)
    xeval = p0[idx]
    count = 0
    while (val - results.min_objective) < stop
        nlow[idx] = xeval
        nhigh[idx] = xeval
        results.opt.lower_bounds = nlow
        results.opt.upper_bounds = nhigh
        p0[idx] = xeval
        init_step = get(kwargs, :init_step, p0 ./ 100.0)
        init_step[idx] = xeval * 1e-6
        init_step[init_step .== 0] .= 0.1
        results.opt.initial_step = init_step
        val, _, _ = optimize!(results.opt, p0)
        likelihood_value = exp(-(val - results.min_objective))
        push!(verb_steps, results.opt.numevals)
        if (likelihood_value == Inf) || (likelihood_value < 0)
            break
        end
        push!(nll_right, likelihood_value)
        push!(x_right, xeval)
        xeval += step
        count += 1
        if count > max_steps
            @warn "Profiling hit limit" max_steps
            break
        end
    end
    verb_count += count
    x = append!(reverse(x_left), x_right[2:end])
    nll = append!(reverse(nll_left), nll_right[2:end])
    if prior !== nothing
        nll = nll .* prior.(x)
    end
    param.likelihood_x = x
    param.likelihood_y = nll
    if verbose
        @debug "Profile steps" verb_count verb_steps
    end
    return x, nll
end

"""
    compute_profiled_uncertainties!(results; kwargs...)

Calculate profiled uncertainties for all parameters in the model.

# Keywords
- `CL::Float64=nothing`: Confidence level
- `σ::Float64=nothing`: Number of sigma
- `mode::String="FC"`: Uncertainty mode ("FC" for Feldman-Cousins)
- `arraystep=nothing`: Per-parameter step sizes
"""
function compute_profiled_uncertainties!(results; kwargs...)
    kwargs = Dict(kwargs)
    CL = get(kwargs, :CL, nothing)
    σ = get(kwargs, :σ, nothing)
    mode = get(kwargs, :mode, "FC")

    arraystep = get(kwargs, :arraystep, nothing)
    if CL !== nothing
        α = CL
    elseif σ !== nothing
        α = erf(σ / sqrt(2))
    else
        α = 0.90
    end
    for (idx, param) in enumerate(results.model.params)
        if arraystep !== nothing
            profile!(param.name, results; step=arraystep[idx], kwargs...)
        else
            profile!(param.name, results; kwargs...)
        end
        uncertainty!(param.name, results; CL=α, mode=mode)
    end
end

"""
    uncertainty!(name, results; CL=nothing, σ=nothing, mode="FC", prior=nothing)

Extract confidence intervals from profiled likelihood for the named parameter.
"""
function uncertainty!(name, results;
                      CL=nothing, σ=nothing, mode="FC", prior=nothing)
    if CL !== nothing
        α = CL
    elseif σ !== nothing
        α = erf(σ / sqrt(2))
    else
        α = erf(1 / sqrt(2)) # 1 Sigma default
    end
    idx = 0
    for (i, p) in enumerate(results.model.params)
        if p.name == name
            idx = i
        end
    end
    param = results.model.params[idx]
    x, y = param.likelihood_x, param.likelihood_y
    if prior !== nothing
        y = y .* prior.(x)
    end
    y = y / sum(y)
    if mode == "FC"
        cy = sortperm(y; rev=true)
        cum_y = cumsum(y[cy])
        cum_x = x[cy]
        region = cum_x[cum_y .<= α]
        left = minimum(region)
        right = maximum(region)
    end
    param.low, param.high = left, right
    return left, right
end
