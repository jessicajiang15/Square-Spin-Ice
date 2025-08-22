using SparseArrays
using Random
using Combinatorics

function kbits(n::Int, k::Int)
    result = String[]
    for bits in combinations(1:n, k)
        s = fill('0', n)
        for bit in bits
            s[bit] = '1'
        end
        push!(result, join(s))
    end
    return result
end

# 
function construct_map_of_half_filling_states(L)
    temp=sort([parse(Int, x; base=2) for x in kbits(L, div(L, 2))], rev=false)
    dict=Dict{Int, Int}()
    count=1
    for t in temp
        dict[t]=count
        count+=1
    end
    return dict
end

# 
function construct_map_of_n_filling_states(L, n)
    temp=sort([parse(Int, x; base=2) for x in kbits(L, n)], rev=false)
    dict=Dict{Int, Int}()
    count=1
    for t in temp
        dict[t]=count
        count+=1
    end
    return dict
end

# Finds all off-diagonal terms that involve hopping from sets of sites "is" and "js"
# is is a list of sites, js is a list of sites, as in, ex. is=[1,2,3...], js=[5, 6...]
# Restrictions: is \cap js = empty, |is|=|js|. i.e. they cannot 
# assumptions: is, js are a valid set of sites and are ordered from lowest to highest index.
# let's say is and js start from index 1
# TESTED: AUGUST 16
function find_possible_pivots_1(is, js, N)
    # list of all the off diagonal coordinates
    index_list=[]
    for state=0:2^N-1
            occupation_i = 0
            occupation_j = 0
            for i in is
                occupation_i+=getTi(i-1, state)
            end
            
            for j in js
                occupation_j+=getTi(j-1, state)
            end
            temp=state
            if (occupation_i==0 && occupation_j==length(js)) || (occupation_i==length(is) && occupation_j==0)
                for i in is
                    temp=flipBit(i-1 ,temp)
                end

                for j in js
                    temp=flipBit(j-1 ,temp)
                end
            end
        if(state < temp)
            push!(index_list, (state+1, temp+1))
        end
            
    end
    return index_list
end


function find_possible_pivots_sector(is, js, the_map)
    # list of all the off diagonal coordinates
    index_list=[]
    states=collect(keys(the_map))
    for state in states
            occupation_i = 0
            occupation_j = 0
            for i in is
                occupation_i+=getTi(i-1, state)
            end
            
            for j in js
                occupation_j+=getTi(j-1, state)
            end
            temp=state
            if (occupation_i==0 && occupation_j==length(js)) || (occupation_i==length(is) && occupation_j==0)
                for i in is
                    temp=flipBit(i-1 ,temp)
                end

                for j in js
                    temp=flipBit(j-1 ,temp)
                end
            end
        if(state < temp)
            push!(index_list, (state+1, temp+1))
        end
            
    end
    return index_list
end

# generates all is, js combinations which involve hopping of n particles at once, for a 1D system of size L, and with periodic boundary
# conditions being toggled by pbc
# returns (is, js)
function generate_all_n_hopping_pairs(n::Int,Nsites::Int)
    result = []
    all_sites = 1:Nsites

    for chosen_sites in combinations(all_sites, 2n)
        # all ways to split chosen_sites into is and js
        for is in combinations(chosen_sites, n)
            js = setdiff(chosen_sites, is)
            # sort both for consistency (optional)
            if(sort(is)<=sort(js))
                push!(result, (sort(is), sort(js)))
            end
        end
    end

    return result
end

function generate_all_n_hopping_pairs_list(n::Int,Nsites::Int)
    result = []
    all_sites = 1:Nsites

    for chosen_sites in combinations(all_sites, 2n)
        # all ways to split chosen_sites into is and js
        for is in combinations(chosen_sites, n)
            js = setdiff(chosen_sites, is)
            # sort both for consistency (optional)
            if(sort(is)<=sort(js))
                push!(result, [sort(is)..., sort(js)...])
            end
        end
    end

    return result
end

# this finds, for n particles hopping, what's the possible indexes in the matrix
# that correspond to each block
# for example, the block corresponding to hopping between 1, 2 will contain the states
# (1, 2), (5,6), and so on.
# this method does so up to hopping between n particles (if n=2, it will generate
# hoppings of 1 particle and hoppings of 2 particles)
function find_all_up_to_order_n_off_diagonal_terms(n, L)
    all_pivots=[]
    for i=1:n
        te=generate_all_n_hopping_pairs(i, L)
        for t in te
            push!(all_pivots,find_possible_pivots_1(t[1],t[2],L))
        end
    end
    return all_pivots
end

# same as above except you also return the corresponding hopping pair
function find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair(n, L)
    all_pivots=[]
    for i=1:n
        te=generate_all_n_hopping_pairs(i, L)
        for t in te
            push!(all_pivots,(t, find_possible_pivots_1(t[1],t[2],L)))
        end
    end
    return all_pivots
end

function find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair_sector(n, L, the_map)
    all_pivots=[]
    for i=1:n
        te=generate_all_n_hopping_pairs(i, L)
        for t in te
            push!(all_pivots,(t, find_possible_pivots_sector(t[1],t[2],the_map)))
        end
    end
    return all_pivots
end




# find all the indicies that correspond to "o" order hopping
# only takes into account of single particle hopping
# running this for o=1...L and gathering all the possible terms corresponds to 1st order as defined in Rademaker's paper
# finds ALL states which allow hoppings between single sites separated by distance o
function find_possible_pivots(o, N, pbc=true)
    # first element: first site hopping to third, second element: second site hopping to fourth
    max_index= (N-o == o) ? N-o : N                  
    index_list=[[] for i=1:max_index]
    for i=0:2^N-1
        for j=1:N
            occupation=getTi(j-1, i)

            # occupation at the neighboring site
            occupation_next=0
            # get its neighbors occupation. If not PBC and last site, there's no neighbor. If PBC, it's neighbors with the first site
            if(j+o<=N)
                occupation_next=getTi(j+o-1, i)
            else
                # if periodic boundary conditions, we get the first site
                if(pbc)
                    occupation_next=getTi(o-(N-j)-1, i)               
                end
            end
            if(occupation_next+occupation==1)
                if(pbc)
                    if(j+o<=N)
                        hopped_state=flipBit(j-1 ,flipBit(j+o-1, i))
                    elseif(N-o != o)
                        hopped_state=flipBit(j-1 ,flipBit(o-(N-j)-1, i))
                    else
                        continue
                    end
                elseif(j+o<=N)
                    hopped_state=flipBit(j-1 ,flipBit(j+o-1, i))
                end
                if(i<hopped_state)
                    push!(index_list[j], (i+1, hopped_state+1))
                end
                
            end
                            
        end
    end
    return index_list
end

# find the next set of pivots (largest)
# order: what order of hopping you have (2- 2nd nearest neighbor hoppings)
# H_int: intermediate Hamiltonian
# uses the maximum possible coefficient
function find_pivot(pivots, H_int)
    maximums=[]
    for piv in pivots
       push!(maximums, maximum(abs.([H_int[p...] for p in piv])))
    end
    return pivots[sortperm(maximums, rev=true)]
end

function find_pivot_first_n(pivots, H_int, n)
    maximums=[]
    for piv in pivots
       push!(maximums, maximum(abs.([H_int[p...] for p in piv])))
    end
    if(n<length(maximums))
        return pivots[sortperm(maximums, rev=true)][1:n]
    else
        return pivots[sortperm(maximums, rev=true)]
    end
end

function find_pivot_first_n_with_maximums(pivots, H_int, n)
    maximums=[]
    for piv in pivots
       push!(maximums, maximum(abs.([H_int[p...] for p in piv])))
    end
    if(n<length(maximums))
        return (maximum(maximums), pivots[sortperm(maximums, rev=true)][1:n])
    else
        return (maximum(maximums), pivots[sortperm(maximums, rev=true)])
    end
end

function find_pivot_first_n_monomials(pivots, H_int, n)
    maximums=[]
    for piv in pivots
       push!(maximums, maximum(abs.(projectors_to_monomials([H_int[p...] for p in piv]))))
    end
    if(n<length(maximums))
        return pivots[sortperm(maximums, rev=true)][1:n]
    else
        return pivots[sortperm(maximums, rev=true)]
    end
end

# in order of the binary numbers
function get_all_hopping_in_one_block(H, max_order, L)
    pivots=find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair(max_order, L)
    temp1=[]
    for piv in pivots
        temp=[]
        for p in piv[2]
            push!(temp, H[p...])
        end
        push!(temp1, (piv[1], temp))
    end
    return temp1
end

function get_all_hopping_in_one_block_monomial(H, max_order, L)
    pivots=find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair(max_order, L)
    temp1=[]
    for piv in pivots
        temp=[]
        for p in piv[2]
            push!(temp, H[p...])
        end
        push!(temp1, (piv[1], projectors_to_monomials(temp)))
    end
    return temp1
end

function get_all_hopping_in_one_block(H, pivots, max_order, L)
    temp1=[]
    for piv in pivots
        temp=[]
        for p in piv[2]
            push!(temp, H[p...])
        end
        push!(temp1, (piv[1], temp))
    end
    return temp1
end

# givens rotation matrix
function givens_rotation_matrix(n,i,j,θ)
    G = Matrix{Float64}(I,(n,n))
    G[i,i] = G[j,j] = cos(θ)
    G[i,j] = sin(θ)
    G[j,i] = -sin(θ)
    return G
end

# number operator at site j, j=1...N
# can rotate back to the original basis to see its structure
function get_number_operator(j, L)
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];
    for i=0:2^L-1
        occupation=getTi(j-1, i)
        if(occupation==1)
            push!(cols, i+1)
            push!(rows, i+1)
            push!(values, 1)
        end
    end
    return sparse(rows, cols, values)
end

function get_number_operator_dense(j, L)
    N=zeros(Float64, 2^L, 2^L)
    for i=0:2^L-1
        occupation=getTi(j-1, i)
        if(occupation==1)
            N[i+1, i+1] = 1
        end
    end
    return N
end

function get_total_number_operator(L)
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];
    for i=0:2^L-1
        push!(cols, i+1)
        push!(rows, i+1)
        push!(values, countBits(n))
    end
    return sparse(rows, cols, values)
end

function get_total_number_operator_dense(L)
    N=zeros(Float64, 2^L, 2^L)
    for i=0:2^L-1
        N[i+1, i+1] = countBits(i)
    end
    return N
end

# includes the site at i_0
function get_N_R_operator(L, i_0)
    N=zeros(Float64, 2^L, 2^L)
    for i=0:2^L-1
        N[i+1, i+1] = countBits(i & (2^i_0 - 1))
    end
    return N
end

function get_N_R_operator(L, map, i_0)
    the_values=collect(keys(map))
    N=zeros(Float64, length(the_values), length(the_values))
    for i in the_values
        index=map[i]
        N[index, index] = countBits(i & (2^i_0 - 1))
    end
    return N
end

function get_N_L_operator(L, map, i_0)
    the_values=collect(keys(map))
    N=zeros(Float64, length(the_values), length(the_values))
    for i in the_values
        index=map[i]
        N[index, index] = countBits(i & (2^L - 2^i_0))
    end
    return N
end

function get_N_R_operator_non_interacting(L, i_0)
    N=zeros(Float64, L, L)
    for i=1:L
        N[i, i] = i<=i_0 ? 1 : 0
    end
    return N
end

function get_N_R_operator_eigenstate_basis(L, eigenstates, i_0)
    N_R=get_N_R_operator(L, i_0)
    N_R=eigenstates'*N_R*eigenstates
    N_R=N_R-Diagonal(diag(N_R))
    return N_R
end

function get_total_N_operator_eigenstate_basis(L, eigenstates)
    N_tot=get_total_number_operator_dense(L)
    N_tot=eigenstates'*N_tot*eigenstates
    return N_tot
end

function states_with_bit1_at(L::Int, i::Int)
    @assert 1 ≤ i ≤ L
    mask = UInt(1) << (i-1)
    lower_max = UInt(1) << (i-1)     # bits below i
    upper_max = UInt(1) << (L-i)     # bits above i
    a=[mask | lower | (upper << i)
        for lower in UInt(0):UInt(lower_max-1),
            upper in UInt(0):UInt(upper_max-1)]
    return Int64.(a)
end

# note that this is index 1
function half_fill_bit1_at(L::Int, i::Int)
    k = div(L, 2)
    @assert 1 ≤ i ≤ L
    @assert k ≥ 1 && k ≤ L "half-filling incompatible with fixed 1 at i"
    maski = UInt(1) << (i-1)
    positions = [p for p in 1:L if p != i]     # remaining positions
    a=[maski | sum(UInt(1) << (p-1) for p in comb)
        for comb in combinations(positions, k-1)]
    return Int64.(a)
end

# in the entire space
function get_n_i_operator(i, L)
    N_i=zeros(Int64, 2^L, 2^L)
    states=states_with_bit1_at(L::Int, i::Int)
    for state in states
        N_i[state+1, state+1] = 1
    end
    return N_i
end

# in the half-filling space
# use one index
function get_n_i_operator_half_filling(i, L, the_map)
    num_states=length(collect(Base.keys(the_map)))
    N_i=zeros(Float64, num_states, num_states)
    possible_states=half_fill_bit1_at(L, i)
    for state in possible_states
        N_i[the_map[state],the_map[state]]=1
    end
    return N_i
end

# get an approximate LIOM by doing a very short time evolution
function get_LIOM_from_n_i(i, H, L, t)
    N_tot=get_n_i_operator(i, L)
    return exp(im * H * t) * N_tot * exp(-im * H * t)
end

# for half-filling
function get_LIOM_from_n_i(i, H, L, t,the_map)
    N_tot=get_n_i_operator_half_filling(i, L,the_map)
    return exp(im * H * t) * N_tot * exp(-im * H * t)
end

function get_N_R_operator_eigenstate_basis(L, the_map,eigenstates, i_0)
    N_R=get_N_R_operator(L,the_map, i_0)
    N_R=eigenstates'*N_R*eigenstates
    N_R=N_R-Diagonal(diag(N_R))
    return N_R
end

function offdiag_norm(A::SparseMatrixCSC; p=2)
    rows, cols, vals = findnz(A)
    offdiag_vals = vals[rows .!= cols]
    return norm(offdiag_vals, p)
end


# assume j starts from index 1
function find_h_j_from_H(H, j)
    return H[2^(j-1)+1,2^(j-1)+1]
end

function hs_to_dict(hs)
    dict=Dict{[Int], Float64}()
    for i=1:length(hs)
        dict[[i]]=hs[i]
    end
    return dict
end

# find all interaction terms up to order (representing the number of particles interacting)
# assume you have hs already
function find_all_interaction_terms(H, L, order, hs)
    # start out with only hs
    previous_Vs = Dict([x] => hs[x] for x=1:length(hs))
    ls=collect(1:L)
    # for each order
    for o=2:order
        # find all of the possible interactions. for example:
        # order 3, possible interactions include 123, 234, ...
        possible_interactions=collect(Combinatorics.combinations(ls, o))
        # each possbile interaction is a set of sites: V_is, for each of these sets of sites is, we need to find V_is
        for is in possible_interactions
            # we find the state corresponding to is (1s at is 0s everywhere else)
            state=find_state_corresponding_to_is(is)

            # get the corresponding diagonal element
            temp=H[state+1, state+1]

            # now for each of the previous orders, we need to subtract from temp
            for a=1:o-1
                # get all possible site combinations at this previous order and subtract them out
                interactions=collect(Combinatorics.combinations(is, a))
                for i in interactions
                    temp-=previous_Vs[i]
                end
            end
            # add it to the dictionary and continue
            previous_Vs[is] = temp
        end
    end
    return previous_Vs
end

# find all interaction terms up to order (representing the number of particles interacting)
# assume you have hs already
function find_all_interaction_terms(H, L, order)
    # start out with only hs
    hs= [find_h_j_from_H(H, i) for i=1:L]
    previous_Vs = Dict([x] => hs[x] for x=1:length(hs))
    ls=collect(1:L)
    # for each order
    for o=2:order
        # find all of the possible interactions. for example:
        # order 3, possible interactions include 123, 234, ...
        possible_interactions=collect(Combinatorics.combinations(ls, o))
        # each possbile interaction is a set of sites: V_is, for each of these sets of sites is, we need to find V_is
        for is in possible_interactions
            # we find the state corresponding to is (1s at is 0s everywhere else)
            state=find_state_corresponding_to_is(is)

            # get the corresponding diagonal element
            temp=H[state+1, state+1]

            # now for each of the previous orders, we need to subtract from temp
            for a=1:o-1
                # get all possible site combinations at this previous order and subtract them out
                interactions=collect(Combinatorics.combinations(is, a))
                for i in interactions
                    temp-=previous_Vs[i]
                end
                
            end
            # add it to the dictionary and continue
            previous_Vs[is] = temp
        end
    end
    return previous_Vs
end

#order < |is|
function find_all_Vi(is, hs)
    # organize this as follows: list of lists. first list, 2 sites. second list, 3 sites, third list, 4 sites etc. 
    # indexing: 1, 2
    previousVs = []
    append!(previousVs, hs)
    for o in orders
        sites=collect(Combinatorics.combinations(is, o))
        for s in sites
            
        end
        find_V_ij_from_H(H, hs, is, previous_V)
    end
end

# assume j starts from index 1
# you need to calculate all of the hs first
# is refers to all of the "interacting sites"
# you need to calculate previous "lower order" Vs
function find_V_ij_from_H(H, is, previous_V)
    state_1=find_state_corresponding_to_is(is)
    temp=H[state_1+1, state_1+1]
    for i in is
        temp = temp - h[i]
    end
    return temp
end

#index 1 for i_s, j_s
function find_state_corresponding_to_is_to_js(i_s, j_s)
    state_1=0
    state_2=0
    for i in i_s
        state_1=flipBit(i-1, state_1)
    end
    for j in j_s
        state_2=flipBit(j-1, state_2)
    end
    return (state_1, state_2)
end

function find_state_corresponding_to_is(i_s)
    state_1=0
    for i in i_s
        state_1=flipBit(i-1, state_1)
    end
    return (state_1)
end

# find the hopping corresponding to between sites i_s and j_s, representing hopping from
# sites i_s to j_s
# must assume that i_s and j_s have no overlap
# also, this is finding the part of t_ij in front of the polynomial of n_i 
# such that you enforce occupations 
function find_t_i_j_from_H(H, is, js)
    state1, state2=find_state_corresponding_to_is_to_js(is, js)
    return abs(H[state1+1, state2+1])
end

# this gives all of the t_ij coeefficients in the i, j block, squared and summed.
function find_norm_t_i_j_from_H(H, is, js)
    state1, state2=find_state_corresponding_to_is_to_js(is, js)
    return abs(H[state1+1, state2+1])
end

function find_all_ts(H, L)
    ts=[]
    for n=1:L-1
        temp=generate_all_n_hopping_pairs(n,L)
        for pair in temp
            push!(ts,find_t_i_j_from_H(H, pair[1], pair[2]))
        end
    end
    return ts
end

# don't add anything less than some tolerance value
function find_all_ts(H, L, tol)
    ts=[]
    for n=1:L-1
        temp=generate_all_n_hopping_pairs(n,L)
        for pair in temp
            off_diagonal_element=find_t_i_j_from_H(H, pair[1], pair[2])
            if(off_diagonal_element > tol)
                push!(ts,off_diagonal_element)
            end
        end
    end
    return ts
end

# maximum out of each block
function find_all_ts_maximum(H, L, pbc=true)
    temp=[]
    all_off_diagonal_indicies=find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair(n, L)
    for pivots in all_off_diagonal_indicies
        distance = pbc ? max_hopping_distance_pbc(pivots[1], L) :  max_hopping_distance_obc(pivots[1])
        list_off_hopping_magnitudes=[H[piv[1],piv[2]] for piv in pivots[2]]
        push!(temp, abs(maximum(abs.(list_off_hopping_magnitudes))))
    end
    return temp
end

function find_all_ts_with_pairs(H, L)
    ts=[]
    for n=1:L-1
        temp=generate_all_n_hopping_pairs(n,L)
        for pair in temp
            push!(ts,(pair, find_t_i_j_from_H(H, pair[1], pair[2])))
        end
    end
    return ts
end

function find_all_h(H, L)
    hs=[]
    for i=1:L
        push!(hs,find_h_j_from_H(H, i))
    end
    return hs
end

function max_pbc_distances(vecs::Vector{Vector{Int}}, L::Int)
    return [length(v) ≥ 2 ?
            maximum(min(abs(i - j), L - abs(i - j)) for (i, j) in combinations(v, 2)) :
            0
            for v in vecs]
end

function max_pbc_distance_0(vecs, L::Int)
    return maximum(min(abs(i - i_0), L - abs(i - i_0)) for i in vecs)
end

function max_obc_distance_from_i0(vecs, i_0)
    return maximum(abs(i - i_0) for i in vecs)
end

function average_obc_distance_from_i0(vecs, i_0)
    return ceil(Int, mean(abs(i - i_0) for i in vecs))
end

function average_pbc_distance_from_i0(vecs, L, i_0)
    return ceil(Int, mean(min(abs(i - i_0), L - abs(i - i_0)) for i in vecs))
end

function range_of_string(sites, L, pbc)
    if(length(sites)==1)
        return 1
    end
    possible_combos=collect(combinations(sites, 2))
    return pbc ? maximum([min(abs(i[1]-i[2]),L-(abs(i[1]-i[2]))) for i in possible_combos])+1 : maximum([abs(i[1]-i[2]) for i in possible_combos])+1
end

# computes the maximum hopping distance between the hopping indicies in A=([pair1],[pair2])
# ex: (1,2,3) and (4,5,6) for L=8 should have maximum hopping distance 4 in PBC, and
# for OBC it should have maximum hopping distance 7
max_hopping_distance_pbc(A, L::Int) =
    maximum( min.(abs.(permutedims(A[1]) .- A[2]), L .- abs.(permutedims(A[1]) .- A[2])) )


function max_hopping_distance_from_i_0_pbc(A, i_0)
    indicies=collect(Iterators.flatten(A))
    return maximum([min.(abs.((a) .- i_0), L .- abs.((a) .- i_0)) for a in indicies])
end


function average_hopping_distance_from_i_0_pbc(A, i_0)
    indicies=collect(Iterators.flatten(A))
    return mean([min.(abs.((a) .- i_0), L .- abs.((a) .- i_0)) for a in indicies])
end
function max_hopping_distance_from_i_0_obc(A, i_0)
    indicies=collect(Iterators.flatten(A))
    return maximum([abs(a-i_0) for a in indicies])
end

function average_hopping_distance_from_i_0_obc(A, i_0)
    indicies=collect(Iterators.flatten(A))
    return mean([abs(a-i_0) for a in indicies])
end

# same as above but for OBC
max_hopping_distance_obc(A) = max(maximum(A[1]) - minimum(A[2]), maximum(A[2]) - minimum(A[1]))

# up to order n (number of particles)
# this is in the basis of n(1-n)(1-n)n... etc
function find_hopping_terms_sorted_by_distance(H, L, n,pbc=true)
    # each index of temp corresponds to a max distance
    temp=[[] for i=1:div(L,2)]
    all_off_diagonal_indicies=find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair(n, L)
    for pivots in all_off_diagonal_indicies
        distance = pbc ? max_hopping_distance_pbc(pivots[1], L) :  max_hopping_distance_obc(pivots[1])
        list_off_hopping_magnitudes=[H[piv[1],piv[2]] for piv in pivots[2]]
        append!(temp[distance], list_off_hopping_magnitudes)
    end
    return temp
end

# finds the "max" hopping term out of a block (which I suppose depends on the "basis" you consider
# for the polynomial)
function find_hopping_norms_sorted_by_distance_maximum(H, L, n,pbc=true)
    # each index of temp corresponds to a max distance
    temp=[[] for i=1:div(L,2)]
    all_off_diagonal_indicies=find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair(n, L)
    for pivots in all_off_diagonal_indicies
        distance = pbc ? max_hopping_distance_pbc(pivots[1], L) :  max_hopping_distance_obc(pivots[1])
        list_off_hopping_magnitudes=[H[piv[1],piv[2]] for piv in pivots[2]]
        #println("hello")
        #println(list_off_hopping_magnitudes)
        append!(temp[distance], abs(maximum(abs.(list_off_hopping_magnitudes))))
    end
    return temp
end

# same as above, but doesnt add any element less than the specified tolerance
function find_hopping_norms_sorted_by_distance_maximum(H, L, n, tol, pbc=true)
    # each index of temp corresponds to a max distance
    temp=[[] for i=1:div(L,2)]
    all_off_diagonal_indicies=find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair(n, L)
    for pivots in all_off_diagonal_indicies
        distance = pbc ? max_hopping_distance_pbc(pivots[1], L) :  max_hopping_distance_obc(pivots[1])
        list_off_hopping_magnitudes=[H[piv[1],piv[2]] for piv in pivots[2] if abs(H[piv[1], piv[2]]) > tol]
        #println(list_off_hopping_magnitudes)
        if(length(list_off_hopping_magnitudes)>0)
            append!(temp[distance], list_off_hopping_magnitudes)
        end
    end
    return temp
end

function find_hopping_norms_sorted_by_distance(H, L, n,pbc=true)
    # each index of temp corresponds to a max distance
    temp=[[] for i=1:div(L,2)]
    all_off_diagonal_indicies=find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair(n, L)
    for pivots in all_off_diagonal_indicies
        distance = pbc ? max_hopping_distance_pbc(pivots[1], L) :  max_hopping_distance_obc(pivots[1])
        list_off_hopping_magnitudes=[H[piv[1],piv[2]] for piv in pivots[2]]
        #println("hello")
        #println(list_off_hopping_magnitudes)
        append!(temp[distance], norm((list_off_hopping_magnitudes)))
    end
    return temp
end

#max_hopping_distance_from_i_0_obc(A, i_0)
function find_hopping_norms_sorted_by_distance_maximum_from_i0(H, L, n, i_0, pbc=false)
    # each index of temp corresponds to a max distance
    temp=[[] for i=1:div(L,2)]
    all_off_diagonal_indicies=find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair(n, L)
    for pivots in all_off_diagonal_indicies
        distance = pbc ? max_hopping_distance_from_i_0_pbc(pivots[1], i_0) :  max_hopping_distance_from_i_0_obc(pivots[1], i_0)
        list_off_hopping_magnitudes=[H[piv[1],piv[2]] for piv in pivots[2]]
        #println("hello")
        #println(list_off_hopping_magnitudes)
        append!(temp[distance], abs(maximum(abs.(list_off_hopping_magnitudes))))
    end
    return temp
end

# same as above but in a specific sector
function find_hopping_norms_sorted_by_distance_maximum_from_i0_sector(H, L, n, i_0, the_map, pbc=false)
    # each index of temp corresponds to a max distance
    num= pbc ? div(L,2) : abs(max(i_0, L-i_0))
    temp=[[] for i=1:num]
    all_off_diagonal_indicies=find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair_sector(n, L, the_map)
    for pivots in all_off_diagonal_indicies
        distance = pbc ? max_hopping_distance_from_i_0_pbc(pivots[1], i_0) :  max_hopping_distance_from_i_0_obc(pivots[1], i_0)
        list_off_hopping_magnitudes=[H[the_map[piv[1]-1],the_map[piv[2]-1]] for piv in pivots[2]]
        #println("hello")
        #println(list_off_hopping_magnitudes)
        append!(temp[distance], abs(maximum(abs.(list_off_hopping_magnitudes))))
    end
    return temp
end

function find_hopping_norms_sorted_by_distance_average_from_i0_sector(H, L, n, i_0, the_map, pbc=false)
    # each index of temp corresponds to a max distance
    num= pbc ? div(L,2) : abs(max(i_0, L-i_0))
    temp=[[] for i=1:num]
    all_off_diagonal_indicies=find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair_sector(n, L, the_map)
    for pivots in all_off_diagonal_indicies
        distance = pbc ? average_hopping_distance_from_i_0_pbc(pivots[1], i_0) :  average_hopping_distance_from_i_0_obc(pivots[1], i_0)
        list_off_hopping_magnitudes=[H[the_map[piv[1]-1],the_map[piv[2]-1]] for piv in pivots[2]]
        append!(temp[ceil(Int, distance)], abs(maximum(abs.(list_off_hopping_magnitudes))))
    end
    return temp
end

function find_hopping_norms_sorted_by_distance_average_from_i0(H, L, n, i_0, pbc=false)
    # each index of temp corresponds to a max distance
    num= pbc ? div(L,2) : abs(max(i_0, L-i_0))
    temp=[[] for i=1:num]
    all_off_diagonal_indicies=find_all_up_to_order_n_off_diagonal_terms_with_hopping_pair(n, L)
    for pivots in all_off_diagonal_indicies
        distance = pbc ? average_hopping_distance_from_i_0_pbc(pivots[1], i_0) :  average_hopping_distance_from_i_0_obc(pivots[1], i_0)
        list_off_hopping_magnitudes=[H[piv[1],piv[2]] for piv in pivots[2]]
        append!(temp[ceil(Int, distance)], abs(maximum(abs.(list_off_hopping_magnitudes))))
    end
    return temp
end

# please note: there is distance zero as well here, whereas the hopping
# function above does not have distance zero
# returns a list where the index is the distance from i_0-1 (average distance, ceiling)
# starting from a distance of zero (so distance zero has index 1)
function find_interaction_terms_sorted_by_distance_average_from_i0(H, L, i_0, pbc=false)
    # each index of temp corresponds to a max distance
    num= pbc ?  abs(max(i_0, L-i_0)) : max(i_0, (L-i_0)+1)
    interactions=[[] for i=1:num]
    all_sites=collect(combinations(1:L))
    for sites in all_sites
        distance=pbc ? average_pbc_distance_from_i0(sites, L, i_0) : average_obc_distance_from_i0(sites, i_0)
        interaction=get_interaction_term(sites, H)
        push!(interactions[distance+1], interaction)
    end
    return interactions
end

# only generate up to order 2 interactions
# lbits are generated in order 1:L (of localization)
# assume you are in the full sector (but can probably restrict to half-filling-sector)
function get_l_bit_model_hamiltonian(lbits, localization_length, pbc, L)
    H=zeros(Float64, 2^L, 2^L)
    d=Normal(1, 1)
    dict = Dict{Vector{Int}, Float64}()
    for i=1:L
        dict[[i]] = rand(d)
        H+=dict[[i]]*lbits[i]
    end
    a=generate_all_n_hopping_pairs_list(1,L)

    for s in a
        site1=s[1]
        site2=s[2]
        distance= pbc ? min(abs(site1-site2), L-abs(site1-site2)) : abs(site1-site2)
        dict[s] = rand(d)*exp(-distance/localization_length)
        H+=dict[s]*lbits[site1]*lbits[site2]
    end
    return (H, dict)
end

function get_l_bit_model_hamiltonian_half_filling(lbits, localization_length, pbc, L, the_map)
    all_keys=Base.keys(the_map)
    H=zeros(Float64, length(all_keys), length(all_keys))
    d=Normal(1, 1)
    dict = Dict{Vector{Int}, Float64}()
    for i=1:L
        dict[[i]] = rand(d)
        H+=dict[[i]]*lbits[i]
    end
    a=generate_all_n_hopping_pairs_list(1,L)

    for s in a
        site1=s[1]
        site2=s[2]
        distance= pbc ? min(abs(site1-site2), L-abs(site1-site2)) : abs(site1-site2)
        dict[s] = rand(d)*exp(-distance/localization_length)
        H+=dict[s]*lbits[site1]*lbits[site2]
    end
    return (H, dict)
end

# sorted by the range of the strings
function find_interaction_terms_sorted_by_range(H, L, pbc=false)
    # each index of temp corresponds to a max distance
    num = pbc ?  div(L,2)+1 : (L)
    interactions=[[] for i=1:num]
    all_sites=collect(combinations(1:L))
    for sites in all_sites
        distance=range_of_string(sites, L, pbc)
        interaction=get_interaction_term(sites, H)
        push!(interactions[distance], interaction)
    end
    return interactions
end

# remove those that are under some tolerance
function find_interaction_terms_sorted_by_range(H, L, tol, pbc=false)
    # each index of temp corresponds to a max distance
    num = pbc ?  div(L,2)+1 : (L)
    interactions=[[] for i=1:num]
    all_sites=collect(combinations(1:L))
    for sites in all_sites
        distance=range_of_string(sites, L, pbc)
        interaction=get_interaction_term(sites, H)
        if(interaction > tol)
            push!(interactions[distance], interaction)
        end
    end
    return interactions
end

# this is starting from distance 0
function find_all_terms_sorted_by_distance_average_from_i_0(H, L, n, i_0, pbc=false)
    all_interactions=find_interaction_terms_sorted_by_distance_average_from_i0(H, L, i_0, pbc)
    all_hoppings=find_hopping_norms_sorted_by_distance_average_from_i0(H, L, n, i_0, pbc)
    count=1
    for hop in all_hoppings
        append!(all_interactions[count+1], hop)
        count+=1
    end
    return all_interactions
end

function get_N_L_operator(L, i_0)
    N=zeros(Float64, 2^L, 2^L)
    for i=0:2^L-1
        N[i+1, i+1] = countBits(i & ((2^i_0 - 1) << i_0 ))
    end
    return N
end

function get_N_L_operator(L, map, i_0)
    the_values=(collect(keys(map)))
    N=zeros(Float64, length(the_values), length(the_values))
    for i in the_values
        index=map[i]
        N[index, index] = countBits(i & ((2^i_0 - 1) << i_0 ))
    end
    return N
end



# THESE TWO FUNCTIONS SWAP COEFFICIENTS IN TWO BASES
# f[mask] are projector (minterm) coeffs, mask in 0:(2^K-1)
# returns a[mask] = monomial coeffs for ∏_{i∈mask} n_i
    function projectors_to_monomials(f::AbstractVector)
        L = length(f)
        @assert L > 0 && L & (L - 1) == 0 "length must be a power of two"
        K = trailing_zeros(L)  # since L == 2^K
        a = float.(f)
        # subset Möbius inversion: for each bit, subtract the value of the subset without that bit
        for b in 0:K-1
            step = 1 << b
            @inbounds for mask in 0:L-1
                if (mask & step) != 0
                    a[mask + 1] -= a[(mask - step) + 1]
                end
            end
        end
        return a
    end
    
    # inverse transform (monomials -> projectors)
    function monomials_to_projectors(a::AbstractVector)
        L = length(a)
        @assert L > 0 && L & (L - 1) == 0 "length must be a power of two"
        K = trailing_zeros(L)
        f = float.(a)
        for b in 0:K-1
            step = 1 << b
            @inbounds for mask in 0:L-1
                if (mask & step) != 0
                    f[mask + 1] += f[(mask - step) + 1]
                end
            end
        end
        return f
    end

    # sites: set of sites that are in the interaction, e.g. sites=[1,2,3] means
    # you are finding the interaction term h_123
    # H: the hamiltonian
    # different way of getting the interactions with FFT 
    function get_interaction_term(sites, H)
        possible_sites=collect(combinations(sites))
        push!(possible_sites,[])
        sum=0
        for s in possible_sites
            state=0
            for i in s
                state = flipBit(i-1,state)
            end
            sum+=H[state+1, state+1]*(-1)^(length(sites)-length(s))
            #println(H[state+1, state+1]*(-1)^(length(sites)-length(s)))
        end
        return sum
    end

    function each_configs_with_ones(L::Int, ones::AbstractVector{Int})
        base = zero(UInt)
        @inbounds for i in ones
            base |= UInt(1) << (i-1)
        end
        free = setdiff(1:L, ones)
        n = length(free)
        return collect(begin
            x = base
            @inbounds for k in 1:n
                if (m >> (k-1)) & 0x1 == 1
                    x |= UInt(1) << (free[k]-1)
                end
            end
            Int(x)
        end for m in UInt(0):(UInt(1)<<n)-UInt(1))
    end