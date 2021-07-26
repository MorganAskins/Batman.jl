# Batman.jl
_A **B**ayesian **A**nalysis **T**oolkit for **M**onitoring
**A**nti-**N**eutrinos, written in [Julia], which provides a statistical
analysis interface geared towards low-count experiments. BATMAN provides many
customizable methods for evaluating confidence intervals and testing
experimental sensitivity._

|                   |                                        |
| :---------------- | :------------------------------------- |
| **Documentation** | [![][stable-img]][stable-url] [![][dev-img]][dev-url] |
| **Development**   | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] |
| **Deployment**    | [![][docker-img]][docker-url] [![][binder-img]][binder-url] |

---

## Installation
With Julia, you can install the Batman.jl package using the built-in
package manager.
```julia-repl
julia> using Pkg
julia> Pkg.add(PackageSpec(url="https://github.com/MorganAskins/Batman.jl"))
```

## Documentation
- [**STABLE**][stable-url] &mdash; **Documentation of the most recently tagged version.**
- [**DEVEL**][dev-url] &mdash; *Documentation of the master branch.*

### Local Documentation
Documentation for Batman is built using [Documenter.jl].  From the `/docs/`
directory, the documentation can be built locally into `/docs/build` via
```bash
julia make.jl
```
this directory is already ignored by git. The file `/docs/build/index.html` can
be opened in a browser to view the newly built documentation.

## Jupyter Notebooks / Lab


---

[binder-img]: https://mybinder.org/badge_logo.svg
[binder-url]: https://mybinder.org/v2/gh/MorganAskins/Batman.jl/master

[travis-img]: https://travis-ci.com/MorganAskins/Batman.jl.svg?branch=master
[travis-url]: https://travis-ci.com/MorganAskins/Batman.jl

[stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[stable-url]: https://MorganAskins.github.io/Batman.jl/stable

[dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[dev-url]: https://MorganAskins.github.io/Batman.jl/dev

[codecov-img]: https://codecov.io/gh/MorganAskins/Batman.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/MorganAskins/Batman.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/sb0jkbj2thwdu3mr/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/MorganAskins/batman-jl/branch/master

[docker-img]: https://img.shields.io/docker/cloud/build/morganaskins/batman.jl
[docker-url]: https://hub.docker.com/r/morganaskins/batman.jl

[Julia]: http://julialang.org/
[Documenter.jl]: https://github.com/JuliaDocs/Documenter.jl
