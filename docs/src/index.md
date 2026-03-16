```@meta
Author = "Morgan Askins"
```
# Batman.jl
The **B**ayesian **A**nalysis **T**oolkit for **M**onitoring
**A**nti-**N**eutrinos, is a package designed to perform statistical
analysis on small datasets in a generic way. The user may decide the
type of analysis to be performed including the posterior and prior
probability distributions, all while naturally including nuisance
parameters (such as systematic uncertainties).

Source code available:
[github.com/MorganAskins/Batman.jl](https://github.com/MorganAskins/Batman.jl)

## Package Features

- Bayesian and Frequentist statistical analysis
- Profile likelihood with Feldman-Cousins confidence intervals
- Poisson counting experiments with systematic uncertainties
- Spectral fitting (1D and multi-dimensional)
- Custom model building with flexible PDF construction
- Kernel density estimation and histogram-based PDFs

## Installation

Batman.jl requires Julia 1.10 or later. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/MorganAskins/Batman.jl")
```

Or using the Pkg REPL:
```
pkg> add https://github.com/MorganAskins/Batman.jl
```

## Quick Start

```julia
using Batman

# Create a counting experiment
m = CountingExperiment()
add_component!(m, "Signal", 20.0)
add_component!(m, "Background", 30.0; σ=5.0)
set_counts!(m, 55)

# Minimize and compute uncertainties
results = minimize!(m)
compute_profiled_uncertainties!(results; σ=1)
pretty_results(results)
```

## Batsignal

Logo for Batman.jl based on:
```julia
using PyPlot
p(x, f) = fill_between(x, f.(x), color="black")

function bat(x; positive=true)
    H(x) = x >= 0 ? 1.0 : 0.0
    σ(x) = @. √(1-x^2.0)
    e(x) = @. 3σ(x/7.0)
    s(x) = @. 4.2 - 0.5x - 2.8σ(0.5x-0.5)
    b(x) = @. σ(abs(2-x)-1)-x.^2/11 + 0.5x - 3
    q(x) = @. 7.5*x*sign(-x) + 8.4
    j(x) = @. -3.0*x*sign(-x) + 0.2
    k(x) = @. 1.7
    if positive
        if abs(x)>7; return 0
        elseif abs(x)>3; return e.(abs.(x))
        elseif abs(x)>1; return s.(abs.(x))
        elseif abs(x)>0.8; return q(x)
        elseif abs(x)>0.5; return j(x)
        else; return k(x)
        end
    else
        if abs(x)>7; return 0
        elseif abs(x)>4; return -e.(abs.(x))
        else; return b(abs.(x))
        end
    end
end

p(-7:0.01:7, a->bat(a; positive=false))
p(-7:0.01:7, a->bat(a; positive=true))
```
![](assets/logo.svg)
