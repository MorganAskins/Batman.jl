abstract type Variable end

mutable struct Parameter <: Variable
  name::String
  info::String
  min::Float64
  max::Float64
  init::Float64
  units::String
  constant::Bool
  bounds::Tuple{Float64, Float64}
  # Post fit information
  samplewidth::Float64
  efficiencies::Dict{Symbol, Float64}
  fit::Float64
  low::Float64
  high::Float64
  likelihood_x::Array{Float64}
  likelihood_y::Array{Float64}
  function Parameter( name::String; kwargs... )
    kwargs = Dict(kwargs)
    info         = get(kwargs, :info, "")
    min          = get(kwargs, :min, -Inf)
    max          = get(kwargs, :max, Inf)
    init         = get(kwargs, :init, 0.0)
    units        = get(kwargs, :units, "")
    constant     = get(kwargs, :constant, false)
    bounds       = get(kwargs, :bounds, (-Inf, Inf) )
    samplewidth  = get(kwargs, :samplewidth, 0.0)
    efficiencies = get(kwargs, :efficiencies, Dict{Symbol, Float64}())
    new( name, info, min, max, init, units, constant, bounds,
         samplewidth, efficiencies, 
         0.0, 0.0, 0.0,
         Array{Float64}(undef, 0),
         Array{Float64}(undef, 0),
       )
  end

end

mutable struct Constant <: Variable
  name::String
  info::String
  init::Float64
  units::String
  constant::Bool
  function Constant( name::String, init::Float64; 
                     info::String="", units::String="" ) 
    new( name, info, init, units, true )
  end
end

struct Dataset <: Variable
  name::String
  init::String
  constant::Bool
  #init::Symbol
  function Dataset( name::String ; info::String="" )
    new( name, name, true )
  end
end

mutable struct Observable{T<:Real}
  name::String
  sname::Symbol
  min::T
  max::T
end
