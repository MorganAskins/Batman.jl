"""
    function profile!(name, results; kwargs...)

Compute the profiled likelihood over a single parameter of
interest given by `name`.

# Arguments
- `name::String`: parameter of interest
- `results::Results`: fit results
- `stop::Float64=5.0`: Maximum ΔLog(Likelihood)
- `step::Float64=1.0`: Incremental step of the parameter of interest
- `prior=nothing`: Prior probability of the parameter of interest
"""
function profile!(name, results; kwargs...)
  kwargs = Dict(kwargs)
  stop      = get(kwargs, :stop, 5.0)
  step      = get(kwargs, :step, 1.0)
  prior     = get(kwargs, :prior, nothing)
  verbose   = get(kwargs, :verbose, false)
  max_steps = get(kwargs, :maxsteps, 1000)
  
  fitterOptions = get(kwargs, :fitteroptions, nothing)
  if fitterOptions != nothing
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
    @warn idx
    return 0, 0
  end
  ## Store these and report if verbose
  verb_count = 0
  verb_steps = Array{Float64}(undef, 0)
  #component = results.model.component_dict[name]
  param = results.model.params[idx]
  ## Go left
  x_left = []
  nll_left = []
  nlow = copy(results.lower_bounds)
  nhigh = copy(results.upper_bounds)
  p0 = copy(results.min_parameters)
  val = copy(results.min_objective)
  xeval = p0[idx]
  count = 0
  while (val-results.min_objective) < stop
    nlow[idx] = xeval
    nhigh[idx] = xeval
    results.opt.lower_bounds = nlow
    results.opt.upper_bounds = nhigh
    p0[idx] = xeval
    init_step = get(kwargs, :init_step, p0./100.0)
    init_step[idx] = xeval * 1e-6
    init_step[(init_step .== 0)] .= 0.1
    results.opt.initial_step=init_step
    val, minx, ret = optimize!(results.opt, p0)
    likelihood_value = exp(-(val-results.min_objective))
    push!(verb_steps, results.opt.numevals)
    if (likelihood_value == Inf) || (likelihood_value < 0)
      println("CYA", likelihood_value)
      break
    end
    push!(nll_left, likelihood_value )
    push!(x_left, xeval)
    xeval -= step
    count += 1
    if count > max_steps 
      @warn "Profiling hit limit" max_steps
      break
    end
  end
  verb_count += count
  ## Go right
  x_right = []
  nll_right = []
  nlow = copy(results.lower_bounds)
  nhigh = copy(results.upper_bounds)
  p0 = copy(results.min_parameters)
  val = copy(results.min_objective)
  xeval = p0[idx]
  count = 0
  while (val-results.min_objective) < stop
    nlow[idx] = xeval
    nhigh[idx] = xeval
    results.opt.lower_bounds = nlow
    results.opt.upper_bounds = nhigh
    p0[idx] = xeval
    init_step = get(kwargs, :init_step, p0./100.0)
    init_step[idx] = xeval * 1e-6
    init_step[(init_step .== 0)] .= 0.1
    results.opt.initial_step=init_step
    val, minx, ret = optimize!(results.opt, p0)
    likelihood_value = exp(-(val-results.min_objective))
    push!(verb_steps, results.opt.numevals)
    if (likelihood_value == Inf) || (likelihood_value < 0)
      break
    end
    push!(nll_right, likelihood_value )
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
  ## update the component with the likelihood results
  if prior != nothing
    nll = nll .* prior.(x)
  end
  param.likelihood_x = x
  param.likelihood_y = nll
  if verbose
    @show verb_count
    @show verb_steps
  end
  return x, nll
end

function compute_profiled_uncertainties!(results; kwargs...)
  kwargs = Dict(kwargs)
  CL = get(kwargs, :CL, nothing)
  σ = get(kwargs, :σ, nothing)
  mode = get(kwargs, :mode, "FC")

  arraystep = get(kwargs, :arraystep, nothing)
  if CL != nothing
    α = CL
  elseif σ != nothing
    α = erf(σ/sqrt(2))
  else
    α = 0.90
  end
  for (idx,param) in enumerate(results.model.params)
    if arraystep != nothing
      profile!(param.name, results; step=arraystep[idx], kwargs...)
    else
      profile!(param.name, results; kwargs...)
    end
    uncertainty!(param.name, results; CL=α, mode=mode)
  end
end

function uncertainty!(name, results; 
                      CL=nothing, σ=nothing, mode="FC", prior=nothing)
  if CL != nothing
    α = CL
  elseif σ != nothing
    α = erf(σ/sqrt(2))
  else
    α = erf(1/sqrt(2)) # 1 Sigma default
  end
  idx = 0
  for (i, p) in enumerate(results.model.params)
    if p.name == name
      idx = i
    end
  end
  param = results.model.params[idx]
  x, y = param.likelihood_x, param.likelihood_y
  # Apply the arbitrary prior function f(x)
  if prior != nothing
    y = y .* prior.(x)
  end
  # Mode can be FC (Feldman-Cousins), Mode-centered, mean-centered, left, right, etc.
  # Return cumulative, x_low, x_high
  y = y / sum(y)
  cumulative = cumsum(y)
  β = 1 - α
  if mode == "FC"
    cy = sortperm(y; rev=true)
    cum_y = cumsum(y[cy])
    cum_x = x[cy]
    region = cum_x[ cum_y .<= α ]
    left = minimum(region)
    right = maximum(region)
  end
  param.low, param.high = left, right
  left, right
end
