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

function get_N_R_operator_eigenstate_basis(L, the_map,eigenstates, i_0)
    N_R=get_N_R_operator(L,the_map, i_0)
    N_R=eigenstates'*N_R*eigenstates
    N_R=N_R-Diagonal(diag(N_R))
    return N_R
end
