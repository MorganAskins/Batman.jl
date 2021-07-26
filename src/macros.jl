"""
@addfunction f

Using `@addfunction` will write in a new function into the
current namespace at read-time.

# Examples
```julia-repl
addfunction lognormalconstraint(x, μ, σ) = (x-μ)^2/2/σ
```
"""
macro addfunction(f)
  :($f)
end

## Is this obsolete?
macro dataset( name, df )
  :( $(esc(name)) = $df )
end

"""
@add_dataset(name, df::DataFrame)

Introduce a new dataset into the namespace at read-time.
"""
function add_dataset(name, df::DataFrame)
  eval(:( $name = $df ) )
end

function add_dataset(df::DataFrame)
  eval(:( $(Symbol(df))=$df ))
end

function add_dataset(name::String, df::DataFrame)
  eval(:( $(Symbol(name)) = $df ))
end

function rename_dataset(name::Symbol, other::Symbol)
  eval(:( $(name) = $(other) ))
end

function add_array(name::Symbol, arr::Array{Float64})
  eval(:( $(name) = $arr ))
end

function add_array(name::String, arr::Array{Float64})
  eval(:( $(Symbol(name)) = $arr ))
end

"""
@add_function(name, function)

Take a user function and pass ownership to a symbol internal
to the Batman module.
"""
function add_function(name, func)
  eval(:( $name = $func ))
end

function add_function(fun)
  eval(:( $(Symbol(fun))=$fun ))
end

function add_function(name::String, func)
  eval(:( $(Symbol(name))=$func ))
end
