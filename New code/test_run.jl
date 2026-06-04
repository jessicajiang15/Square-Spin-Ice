include("unitTests.jl")

W = parse(Int, ARGS[1])
JJ = parse(Float64, ARGS[2])

println("W = $W")
println("J = $JJ")

ctc_calculation_3(W, JJ)