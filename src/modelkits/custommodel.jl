mutable struct CustomModel <: Model
    params::Vector
    pdflist::Vector{NLogPDF}
    nll::NLogLikelihood
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    function CustomModel()
        new([], NLogPDF[], NLogLikelihood(), Float64[], Float64[])
    end
end

function add_parameters!(m::CustomModel, params)
    for p in params
        push!(m.params, p)
    end
end

function add_nlogpdfs!(m::CustomModel, pdflist::Vector{NLogPDF})
    for pdf in pdflist
        push!(m.pdflist, pdf)
    end
end

function minimize!(m::CustomModel; options=Dict())
    likelihood = NLogLikelihood(m.pdflist)
    add_likelihood!(m, likelihood)
    optimize_model!(m, likelihood; options=options)
end
