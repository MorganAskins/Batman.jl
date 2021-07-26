function rowselect(df::DataFrame, s::Expr)
  opera = String(s.args[1])
  if !occursin(".", opera)
    opera = "."*opera
  end
  opera = Symbol(opera)
  obs = s.args[2]
  val = s.args[3]
  df[eval(Expr(:call, opera, df[!, obs], val)), :]
end
