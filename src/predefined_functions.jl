"""
    lognormal(x, μ, σ)

Compute the negative log of a normal distribution (up to a constant):

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
    (x - μ)^2 / 2 / σ^2
end

"""
    logpoisson(n, x...)

Compute the negative log of a Poisson distribution (up to a constant):

```math
λ = \\sum(x) \\\\
λ - n \\log(λ) + n \\log(n) - n
```
"""
function logpoisson(n, x...)
    λ = sum(x)
    if λ <= 0
        return Inf
    end
    λ - n * log(λ) + n * log(n) - n
end

"""
    batlog(x; up=1000.0, down=-1000.0)

Compute the natural log with bounds to handle very large numbers and zeros.

Returns `down` if `x <= 0`, and clamps the result between `down` and `up`.
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
    clamp(y, down, up)
end
