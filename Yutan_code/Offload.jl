using LinearAlgebra
using SparseArrays
using Arpack
using JLD
using Dates
using KrylovKit

simulation = "E:/UC Davis/Research/Square spin ice/Square-Spin-Ice/Yutan_code/Results/Hamiltonian_even_odd_L=4__J=1.0__h=1.0_Date_Thu_17_Jun_2021_22_39_55.jld"

H_even = load(simulation)["sim"][1]
H_odd = load(simulation)["sim"][2]

println(eigsolve(H_even, 200, :SR, krylovdim=200)[1])
