
include("renormalizationGroupHelpers.jl")
using Convex, LinearAlgebra, SCS
# assume the eigenstates are in the computational basis for the half-filling sector
# prepares the state with all particles on the right in the eigenstate basis
function prepare_half_filling_state(eigenstates, map_states, L, i_0)
    # the binary number corresponding to 111...00
    half_filling_binary_no=2^(i_0)-1
    
    # the half filling state but rewritten in the eigenstate basis
    half_filling_state=[]
    index=map_states[half_filling_binary_no]
    for i=1:length(eigenstates[:, 1])
        push!(half_filling_state,eigenstates[index,i])
    end
    return half_filling_state
end

# prepares a state 000...111
function prepare_half_filling_state(eigenstates, map_states, L, i_0)
    # the binary number corresponding to 111...00
    half_filling_binary_no=2^(i_0)-1

    # index 1    
    # the half filling state but rewritten in the eigenstate basis
    half_filling_state=[]
    index=map_states[half_filling_binary_no]
    for i=1:length(eigenstates[:, 1])
        push!(half_filling_state,eigenstates[index,i])
    end
    return half_filling_state
end

# assumes "state" is written in the eigenbasis
# map_states is a map from state -> index. Fine as long as its consistent
function measure_number_particles_right_of_i_0(state, the_map, eigenstates, L, i_0)
    N_R=get_N_R_operator(L,the_map, i_0)
    N_R=eigenstates'*N_R*eigenstates
    return state' * N_R * state
end

function measure_number_particles_left_of_i_0(state, the_map, eigenstates, L, i_0)
    N_R=get_N_L_operator(L,the_map, i_0)
    N_R=eigenstates'*N_R*eigenstates
    return state' * N_R * state
end

# state is written in the computational basis
# we want to rewrite it in the eigenstate basis
# eigenstates is a column wise list of eigenvectors written in the computational basis
function prepare_in_eigenstate_basis(state, eigenstates)
    return eigenstates' * state
end

function optimize_nr(M)
    A = copy(M)
    n = size(A, 1)
    N_R=get_N_L_operator(L,the_map, div(L,2))
    N_R=eigenvectors'*N_R*eigenvectors
    
    # Decision variables
    d = Variable(size(A, 1)) 
    α = Variable()          # lower eigenvalue bound
    β = Variable()          # upper eigenvalue bound
    t = Variable()          # spread (β - α)
    
    D = diagm(d)
    
    constraints = [
        A + D ⪯ β*I(n),     # largest eigenvalue ≤ β
        α*I(n) ⪯ A + D,     # smallest eigenvalue ≥ α
        t == β - α
    ]
    
    problem = minimize(t, constraints)
    solve!(problem, SCS.Optimizer; silent=true)
    
    # Construct the diagonal matrix D
    D_opt = diagm(Convex.evaluate(d))
    return Convex.evaluate(t)
end