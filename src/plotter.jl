# Plotting functions
using Random

function interval_plot(results, name)
  #comp = results.model.component_dict[name]
  ### TODO fix this
  idx = 0
  for (i, p) in enumerate(results.model.params)
    if p.name == name
      idx = i
    end
  end
  comp = results.model.params[idx]
  # Check if these are already computed, if not compute
  if length(comp.likelihood_x) == 0
    compute_profiled_uncertainties!(results, 0.68)
  end
  x, y = comp.likelihood_x, comp.likelihood_y
  low, high = comp.low, comp.high
  mid = comp.fit
  ns, ps = abs(mid - low), abs(mid - high)
  labeluse = "$name: " *
    " \$\\mu = $(floor(mid*100)/100)" *
    "^{+$(floor(ps*100)/100)}" *
    "_{-$(ns)}\$"
  plot(x, y, label=labeluse)
  vlines(low, 0, y[x .>= low][1], color="red")
  vlines(high, 0, y[x .>= high][1], color="red")
  ylim(0, 1.2)
  xlim(minimum(x), maximum(x))
  legend(loc=2, frameon=false)
end

function _mean_comp(x, y)
  sum(x.*y)
end

function _variance(x, y, μ)
  sum( (x .- μ).^2 .*y )
end

function _covariance(x, y, μx, μy, z)
  N = sum(z)
  k = transpose(x .- μx).*(y .- μy).*z
  sum( k )/N
end

function correlation_plots(results; kwargs...)
  kwargs = Dict(kwargs)
  steps = get(kwargs, :steps, 200)
  figwidth = get(kwargs, :figwidth, 10)
  ## Everything stored in results
  M = length(results.model.params)
  fig, ax = subplots(nrows=M, ncols=M, figsize=(M*figwidth, M*figwidth))

  ## We will be lazy for now and robust later. Lazy way, assume we know the range
  # What does a slice in space look like?
  # We can optimize again, with a constraint in place (bounds)
  levels = exp.(-[Inf,11.83,6.18,2.3,0].^0.5./2)

  for (i, v) in enumerate(results.model.params)
    for (j, w) in enumerate(results.model.params)
      if j < M
        ax[j, i].set_xticklabels([])
      end
      if i > 1
        ax[j, i].set_yticklabels([])
      end
      if i > j
        ax[j, i].axis("off")
        continue
      end
      # profile this parameter
      nlow = copy(results.lower_bounds)
      nhigh = copy(results.upper_bounds)
      p0 = copy(results.min_parameters)
      nll = []
      if i != j
        A = range(minimum(v.likelihood_x), length=steps, stop=maximum(v.likelihood_x)) |> collect
        B = range(minimum(w.likelihood_x), length=steps+5, stop=maximum(w.likelihood_x)) |> collect
        function getnll(a, b)
          nlow[i], nlow[j] = a, b
          nhigh[i], nhigh[j] = a, b
          results.opt.lower_bounds = nlow
          results.opt.upper_bounds = nhigh
          #p0 = copy(results.min_parameters) << This is slow
          p0[i], p0[j] = a, b
          minf, minx, ret = optimize!(results.opt, p0)
          if minf == 0
            minf = 10000
          end
          exp(-minf)
        end
        Z = [getnll(a, b) for b in B, a in A]
        Z = Z / maximum(Z)

        ## Given Z, A, B
        yb = sum(Z, dims=2)
        ya = transpose( sum(Z, dims=1) )
        ya = ya / sum(ya)
        yb = yb / sum(yb)

        mean_a = _mean_comp(A, ya)
        mean_b = _mean_comp(B, yb)
        var_a = sqrt( _variance(A, ya, mean_a) )
        var_b = sqrt( _variance(B, yb, mean_b) )
        cov = _covariance(A, B, mean_a, mean_b, Z)
        cor = cov / var_a / var_b
        ax[i, j].text(0,0.5, @sprintf("Cor: %0.3f", cor))
        ax[j,i].contourf(A, B, Z, levels=levels)
        continue
      end
      ax[i, i].plot( v.likelihood_x, v.likelihood_y )
      ax[i, i].set_title(v.name)
    ax[i,i].set_ylim(0, 1.2)
    end
  end
end

function correlation_plots2(results; kwargs...)
  kwargs = Dict(kwargs)
  steps = get(kwargs, :steps, 200)
  figwidth = get(kwargs, :figwidth, 10)
  ## Everything stored in results
  M = length(results.model.params)
  fig, ax = subplots(nrows=M, ncols=M, figsize=(M*figwidth, M*figwidth))

  levels = exp.(-[Inf,11.83,6.18,2.3,0].^0.5./2)

  for (i, v) in enumerate(results.model.params)
    for (j, w) in enumerate(results.model.params)
      if j < M
        ax[j, i].set_xticklabels([])
      end
      if i > 1
        ax[j, i].set_yticklabels([])
      end
      if i > j
        ax[j, i].axis("off")
        continue
      end
      # profile this parameter
      nlow = copy(results.lower_bounds)
      nhigh = copy(results.upper_bounds)
      p0 = copy(results.min_parameters)
      nll = []
      if i != j
        ## Here we need to go quadrant by quadrant
        xl = range( v.fit, length=steps, stop=minimum(v.likelihood_x) ) |> collect
        xr = range( v.fit, length=steps, stop=maximum(v.likelihood_x) ) |> collect
        yu = range( w.fit, length=steps, stop=maximum(w.likelihood_x) ) |> collect
        yd = range( w.fit, length=steps, stop=minimum(w.likelihood_x) ) |> collect

        function getz(xset, yset)
          # We always start at the best fit, so lets reset p0 to best fit
          coords = []
          p0 = copy(results.min_parameters)
          minf0 = copy(results.min_objective)
          pHold = copy(p0)
          for y in yset
            rowz = Array{Float64}(undef, 0)
            p0 = copy(pHold)
            reset = true
            for x in xset
              nlow[i], nlow[j] = x, y
              nhigh[i], nhigh[j] = x, y
              results.opt.lower_bounds = nlow
              results.opt.upper_bounds = nhigh
              p0[i], p0[j] = x, y
              minf, minx, ret = optimize!(results.opt, p0)
              if minf == 0
                minf = 10000
              end
              push!(rowz, exp(-(minf - minf0)) )
              if reset
                reset = false
                pHold = copy(p0)
              end
            end
            push!(coords, rowz)
          end
          return coords
        end
        z1 = getz(xr, yu)
        z2 = getz(xl, yu)
        z3 = getz(xl, yd)
        z4 = getz(xr, yd)
        # Reorder
        z2 = [reverse(a) for a in z2]
        z3 = reverse([reverse(a) for a in z3])
        z4 = reverse(z4)
        zleft = vcat(z3, z2)
        zright = vcat(z4, z1)
        X = vcat(reverse(xl), xr)
        Y = vcat(reverse(yd), yu)
        Z = Array([vcat(i, j) for (i,j) in zip(zleft, zright)])
        Z = hcat(Z...)'
        Z = Z / maximum(Z)

        ## Given Z, A, B
        yb = sum(Z, dims=2)
        ya = transpose( sum(Z, dims=1) )
        ya = ya / sum(ya)
        yb = yb / sum(yb)

        mean_a = _mean_comp(X, ya)
        mean_b = _mean_comp(Y, yb)
        var_a = sqrt( _variance(X, ya, mean_a) )
        var_b = sqrt( _variance(Y, yb, mean_b) )
        cov = _covariance(X, Y, mean_a, mean_b, Z)
        cor = cov / var_a / var_b
        ax[i, j].text(0,0.5, @sprintf("Cor: %0.3f", cor))
        ax[j,i].contourf(X, Y, Z, levels=levels)
        continue
      end
      ax[i, i].plot( v.likelihood_x, v.likelihood_y )
      ax[i, i].set_title(v.name)
    ax[i,i].set_ylim(0, 1.2)
    end
  end
end
