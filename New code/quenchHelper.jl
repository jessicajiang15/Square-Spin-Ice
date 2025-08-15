
include("renormalizationGroupHelpers.jl")
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