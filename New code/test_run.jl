include("unitTests.jl")

W = parse(Int, ARGS[1])
JJ = parse(Int, ARGS[2])

println("W = $W")
println("J = $JJ")

ctc_calculation_3(W, JJ)