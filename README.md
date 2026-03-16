# Batman.jl

_A **B**ayesian **A**nalysis **T**oolkit for **M**onitoring
**A**nti-**N**eutrinos, written in [Julia], which provides a statistical
analysis interface geared towards low-count experiments. BATMAN provides many
customizable methods for evaluating confidence intervals and testing
experimental sensitivity._

|                   |                                        |
| :---------------- | :------------------------------------- |
| **Documentation** | [![][stable-img]][stable-url] [![][dev-img]][dev-url] |
| **Build Status**  | [![][ci-img]][ci-url] [![][codecov-img]][codecov-url] |

---

## Requirements

- Julia 1.10 or later

## Installation

With Julia, you can install the Batman.jl package using the built-in
package manager.
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

# Create a Poisson counting experiment
m = CountingExperiment()
add_component!(m, "Signal", 20.0)
add_component!(m, "Background", 30.0; σ=5.0)
set_counts!(m, 55)

# Fit and compute uncertainties
results = minimize!(m)
compute_profiled_uncertainties!(results; σ=1)
pretty_results(results)
```

## Documentation

- [**STABLE**][stable-url] &mdash; **Documentation of the most recently tagged version.**
- [**DEV**][dev-url] &mdash; *Documentation of the development branch.*

### Local Documentation

Documentation for Batman is built using [Documenter.jl]. From the `/docs/`
directory, the documentation can be built locally into `/docs/build` via
```bash
julia --project=docs/ docs/make.jl
```
This directory is already ignored by git. Open `/docs/build/index.html` in a
browser to view the newly built documentation.

## License

MIT License - see [LICENSE.md](LICENSE.md) for details.

---

[ci-img]: https://github.com/MorganAskins/Batman.jl/actions/workflows/CI.yml/badge.svg?branch=master
[ci-url]: https://github.com/MorganAskins/Batman.jl/actions/workflows/CI.yml?query=branch%3Amaster

[stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[stable-url]: https://MorganAskins.github.io/Batman.jl/stable

[dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[dev-url]: https://MorganAskins.github.io/Batman.jl/dev

[codecov-img]: https://codecov.io/gh/MorganAskins/Batman.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/MorganAskins/Batman.jl

[Julia]: https://julialang.org/
[Documenter.jl]: https://github.com/JuliaDocs/Documenter.jl
