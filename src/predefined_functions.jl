"""
    lognormal(x, μ, σ)

```math
\\frac{(x-μ)^2}{2σ^2}
```
"""
function lognormal(x, μ, σ)
  if σ == Inf
    return 0
  end
  if σ ≈ 0
    return Inf
  end
  (x-μ)^2/2/σ^2
end

"""
    logpoisson(n, x...)

```math
  λ > 0 \\
  λ = sum([x...])
  λ - n*log(λ) + n*log(n) - n
```math
"""
function logpoisson(n, x...)
  λ = sum([x...])
  if λ <= 0
    return Inf
  end
  λ - n*log(λ) + n*log(n) - n
end

"""
    batlog(x)

Compute the natural log, but create bounds to ignore very large
numbers / infinities.
```jldoctest
julia> batlog.([12.0, 1e6, 0])
3-element Array{Float64, 1}:
  2.4849
  13.8155
  -100000.0
```
"""
function batlog(x::Float64; up::Float64=1000.0, down::Float64=-1000.0)
  if x <= 0
    return -1e8
  end
  y = log(x)
  if y < down
    @debug "Batlog hit floor" y
  end
  if y > up
    @debug "Batlog hit ceiling" y
  end
  minimum([up, maximum([down, y])])
end
