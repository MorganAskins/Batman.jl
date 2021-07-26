mutable struct CustomModel <: Model
  params
  pdflist::Array{NLogPDF}
  nll::NLogLikelihood
  lower_bounds::Array{Float64}
  upper_bounds::Array{Float64}
  function CustomModel()
    new( [], 
        Array{NLogPDF}(undef, 0),
        NLogLikelihood(),
        Array{Float64}(undef, 0),
        Array{Float64}(undef, 0)
       )
  end
end

function add_parameters!( m::CustomModel, params )
  for (idx,p) in enumerate(params)
    push!(m.params, p)
  end
end

function add_nlogpdfs!( m::CustomModel, pdflist::Array{NLogPDF} )
  for pdf in pdflist
    push!(m.pdflist, pdf)
  end
end

function minimize!( m::CustomModel; options=Dict() )
  likelihood = NLogLikelihood( m.pdflist )
  add_likelihood!( m, likelihood )
  optimize_model!( m, likelihood; options=options )
end
