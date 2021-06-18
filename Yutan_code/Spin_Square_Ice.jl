using LinearAlgebra
using SparseArrays
using Arpack
using JLD
using Dates

function coordinate(n;L::Int64=2)
    num_sites = L^2
    i::Int64 = Int(ceil(n/L))
    j::Int64 = mod1(n,L)  #site (i,j) is at i-th row, j-th column
    return (i,j)
end

function bit_pos(coordinate::Tuple{Int64,Int64};L::Int64=2)
    n = (coordinate[1]-1)*L + coordinate[2]
    return n
end

function neib(n::Int64;L::Int64=2)
    coord = coordinate(n,L=L)
    neibs = Tuple{Int64,Int64}[]
    push!(neibs, (mod1(coord[1]+1,L), coord[2]))
    push!(neibs, (mod1(coord[1]-1,L), coord[2]))
    push!(neibs, (coord[1], mod1(coord[2]+1,L)))
    push!(neibs, (coord[1], mod1(coord[2]-1,L)))
    if iseven(coord[1]+coord[2])
        push!(neibs, (mod1(coord[1]+1,L), mod1(coord[2]-1,L)))
        push!(neibs, (mod1(coord[1]-1,L), mod1(coord[2]+1,L)))
    else
        push!(neibs, (mod1(coord[1]+1,L), mod1(coord[2]+1,L)))
        push!(neibs, (mod1(coord[1]-1,L), mod1(coord[2]-1,L)))
    end
    #=convert coordinations to positions in bits=#
    neibs_bit_pos = Set{Int64}()
    for neib in neibs
        push!(neibs_bit_pos, bit_pos(neib,L=L))
    end
    return neibs_bit_pos
end

function neib_list_gen(;L::Int64=L)
    neib_list = Set{Int64}[]
    for n in 1:L^2
        push!(neib_list, neib(n, L=L))
    end
    return neib_list
end

function chk_parity(state::Int64)
    state_binary = digits!(zeros(Int64, 64), state, base = 2)
    if iseven(sum(state_binary))
        return :even
    else
        return :odd
    end
end

function parity_states_list(;L::Int64=L)
    even_state = Int64[]
    odd_state = Int64[]
    even_state_num = Dict{Int64, Int64}()
    odd_state_num = Dict{Int64, Int64}()
    even_state_tot = 0
    odd_state_tot = 0
    for state in 0:(2^(L^2)-1)
        if chk_parity(state) == :even
            even_state_tot += 1
            push!(even_state, state)
            even_state_num[state] = even_state_tot
        else
            odd_state_tot += 1
            push!(odd_state, state)
            odd_state_num[state] = odd_state_tot
        end
    end
    return Dict{Symbol, Any}(:even_state => even_state, :odd_state => odd_state, :even_state_num => even_state_num, :odd_state_num => odd_state_num, :even_state_tot => even_state_tot, :odd_state_tot => odd_state_tot)
end

function Hamiltonian(;L::Int64=2, J=1, h=1, neib_list, state_list, state_num, state_tot)
    H = spzeros(state_tot,state_tot)
    for state in state_list #loop over all states
        state_binary = digits!(zeros(Int64, 64), state, base = 2)
        for i in 1:L^2 #loop over all sites i in a given state
            #diagonal terms interact with h
            if state_binary[i] == 1
                H[state_num[state],state_num[state]] += h/2
            else
                H[state_num[state],state_num[state]] -= h/2
            end
            # off diagonal terms come from flipping bonds, j is neighbor of i
            for j in neib_list[i]
                fliped_state = state ⊻ (1<<(i-1))
                fliped_state = fliped_state ⊻ (1<<(j-1))
                H[state_num[state],state_num[fliped_state]] = J/4
            end
        end
    end
    return H
end

function H_driver(;L::Int64=L, J::Float64=J, h::Float64=h)
    neib_list = neib_list_gen(L=L)
    state_gen = parity_states_list(L=L)
    even_state = state_gen[:even_state]
    even_state_num = state_gen[:even_state_num]
    odd_state_num  = state_gen[:odd_state_num]
    odd_state = state_gen[:odd_state]
    even_state_tot = state_gen[:even_state_tot]
    odd_state_tot = state_gen[:odd_state_tot]
    #now constructine Hamiltonian of the two sectors
    Hamiltonian_even = Hamiltonian(;L=L, J=J, h=h, neib_list=neib_list, state_list=even_state, state_num=even_state_num, state_tot=even_state_tot)
    Hamiltonian_odd = Hamiltonian(;L=L, J=J, h=h, neib_list=neib_list, state_list=odd_state, state_num=odd_state_num, state_tot=odd_state_tot)
    return Dict(:even => Hamiltonian_even, :odd => Hamiltonian_odd)
end

L=4;J=1.0;h=1.0
@time H_blocks = H_driver(L=L, J=J, h=h)

my_time = Dates.now()

time_finished = "Date_$(Dates.format(my_time, "e_dd_u_yyyy_HH_MM_SS"))"
content = "Hamiltonian_even_odd"
save_path = "/nfs/home/zyt329/Research/Square_spin_ice/result/"
save_name = save_path*content*"_L=$(L)__J=$(J)__h=$(h)_"*time_finished*".jld"

save(save_name, "sim", [H_blocks[:even],H_blocks[:odd]])
