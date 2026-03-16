function Base.show(io::IO, m::Results)
    compact = get(io, :compact, false)
    if !compact
        println(io, "Fit converged after $(m.iterations) iterations")
        println(io, "Best fit at: -LnL = $(m.min_objective)")
        println(io, "Values of $(m.min_parameters)")
    end
end

function pretty_results(r::Results)
    df = DataFrame(
        Name=String[],
        Fit=Float64[],
        Interval_Low=Float64[],
        Interval_High=Float64[],
    )
    for v in r.model.params
        push!(df, (v.name, v.fit, v.low, v.high))
    end
    return df
end
