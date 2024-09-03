using TensorOperations
using NearestNeighbors
using Distances
using StaticArrays
using LinearAlgebra
using KrylovKit
using Arpack

function calculate_Vijkl(i, j, k, l, eigenvectors, bonds)
    sum=0
    for b in bonds
        bond1::Int=b.site1.num;
        bond2::Int=b.site2.num;
        sum+=(eigenvectors[:,i][bond1]*eigenvectors[:,j][bond2]*conj(eigenvectors[:,k][bond1])*conj(eigenvectors[:,l][bond2]))
        sum-=(eigenvectors[:,i][bond2]*eigenvectors[:,j][bond1]*conj(eigenvectors[:,k][bond1])*conj(eigenvectors[:,l][bond2]))
        sum-=(eigenvectors[:,i][bond1]*eigenvectors[:,j][bond2]*conj(eigenvectors[:,k][bond2])*conj(eigenvectors[:,l][bond1]))
        sum+=(eigenvectors[:,i][bond2]*eigenvectors[:,j][bond1]*conj(eigenvectors[:,k][bond2])*conj(eigenvectors[:,l][bond1]))
    end
    return -sum/4
end

#this assumes bonds are consecutive sites in 1D 
function build_Vijkl_matrix(eigenvectors)
    v_rolled=circshift(eigenvectors, (1, 0))
    n=size(eigenvectors, 1)
    identity = zeros(Float64, n, n, n)
    for i in 1:size(eigenvectors, 1)
        identity[i, i, i] = 1.0
    end
    V = @tensor W[a, b, c, d] := eigenvectors[i, a] * identity[i,j,k] * (v_rolled)[j, b] *conj.((eigenvectors))[l, c] * identity[l,m,k] *  conj.((v_rolled))[m, d]
    V1 = permutedims(V, (2, 1, 3, 4))
    V2 = permutedims(V, (1, 2, 4, 3))
    V3 = permutedims(V, (2, 1, 4, 3))
    V=V-V1-V2+V3
    #println(V)
    temp=ifelse.(0.25 * abs.(V) .< 1e-12, 0, -0.25*V)
    return temp
end

function build_Vijkl_matrix_no_symmetry(eigenvectors, eigenvalues, bonds)
    temp=zeros(Float64,length(eigenvalues),length(eigenvalues),length(eigenvalues),length(eigenvalues))
    for i=1:length(eigenvalues)
        for j=1:length(eigenvalues)
            for k=1:length(eigenvalues)
                for l=1:length(eigenvalues)
                    temp[i,j,k,l]=calculate_Vijkl(i, j, k, l, eigenvectors, bonds)
                end
            end
        end
        return temp
    end
    return ifelse.(abs.(V) .< 1e-12, 0, V)
end

function build_energy_differences_matrix(eigenvalues)
    temp=zeros(Float64,length(eigenvalues),length(eigenvalues),length(eigenvalues),length(eigenvalues))
    Threads.@threads for i=1:length(eigenvalues)
        Threads.@threads for j=1:i
            te=eigenvalues[i]+eigenvalues[j]
            Threads.@threads for k=1:length(eigenvalues)
                Threads.@threads for l=1:k
                    value=te-eigenvalues[k]-eigenvalues[l]
                    temp[i,j,k,l]=value
                    temp[i,j,l,k]=value
                    temp[j,i,l,k]=value
                    temp[j,i,k,l]=value


                    temp[k,l,i,j]=-value
                    temp[l,k,i,j]=-value
                    temp[l,k,j,i]=-value
                    temp[k,l,j,i]=-value
                end
            end
        end
    end
    return ifelse.(abs.(temp) .< 1e-12, 1.0, temp)
end


function build_energy_differences_matrix_no_symmetry(eigenvalues)
    temp=zeros(Float64,length(eigenvalues),length(eigenvalues),length(eigenvalues),length(eigenvalues))
    Threads.@threads for i=1:length(eigenvalues)
        Threads.@threads for j=1:length(eigenvalues)
            te=eigenvalues[i]+eigenvalues[j]
            Threads.@threads for k=1:length(eigenvalues)
                Threads.@threads for l=1:length(eigenvalues)
                    value=te-eigenvalues[k]-eigenvalues[l]
                    temp[i,j,k,l]=value
                end
            end
        end
    end
    return ifelse.(abs.(temp) .< 1e-12, 1, temp)
end

function build_W_matrix_tensor_op(eigenvalues, eigenvectors)
    V=build_Vijkl_matrix(eigenvectors)
    energy=build_energy_differences_matrix(eigenvalues)
    return V./energy
end

function delta(i,j)
    return i==j ? 1 : 0
end

function calculate_Wijkl(i,j,k,l,eigenvectors,eigenvalues, bonds)
    vijkl=calculate_Vijkl(i, j, k, l, eigenvectors, bonds)
    if(i==j==l==k||(abs(eigenvalues[i]+eigenvalues[j]-eigenvalues[k]-eigenvalues[l])<1e-20)||i==k&&j==l||i==l&&j==k)
        return 0
    end
    return vijkl/(eigenvalues[i]+eigenvalues[j]-eigenvalues[k]-eigenvalues[l])
end

function construct_W_matrix(eigenvectors,eigenvalues, bonds)
    W=zeros(Float64,length(eigenvalues),length(eigenvalues),length(eigenvalues),length(eigenvalues))
    for i=1:length(eigenvalues)
        for j=1:length(eigenvalues)
            for k=1:length(eigenvalues)
                for l=1:length(eigenvalues)
                    W[i,j,k,l]=calculate_Wijkl(i,j,k,l,eigenvectors,eigenvalues, bonds)
                end
            end
        end
    end
    return W
end

function construct_W_matrix_symmetry(eigenvectors,eigenvalues, bonds)
    W=zeros(Float64,length(eigenvalues),length(eigenvalues),length(eigenvalues),length(eigenvalues))
    Threads.@threads for j=1:length(eigenvalues)
        Threads.@threads for i=1:j
            Threads.@threads for l=1:length(eigenvalues)
                Threads.@threads for k=1:l
                    W[i,j,k,l]=calculate_Wijkl(i,j,k,l,eigenvectors,eigenvalues, bonds)
                    W[j, i, l, k]=W[i,j,k,l]
                    W[l,k,j, i]=conj(W[i,j,k,l])
                    W[k, l, i, j]=conj(W[i,j,k,l])
                    W[i, j, l, k]=-(W[i,j,k,l])
                    W[j, i, k, l]=-W[i,j,k,l]
                    W[l, k, j, i]=-conj(W[i,j,k,l])
                    W[k, l, i, j]=-conj(W[i,j,k,l])
                end
            end
        end
    end
    return W
end

#assume operator is in matrix form
function find_norm_operator(W)
    term1=0
    for i=1:(size(W, 1))
        for j=1:(size(W, 1))
            term1+=W[i,j,i,j]
        end
    end
    term2=0
    for i=1:(size(W, 1))
        for j=1:(size(W, 1))
            temp=0
            for k=1:(size(W, 1))
                temp+=abs(W[k,i,k,j])
            end
            term2+=abs2(temp)
        end
    end
    term3=norm(W)
    return (0.25)*term1^2+term2+0.25*term3
end

function find_norm_operator(W)
    n = size(W, 1)
    
    term1 = sum(W[i, j, i, j] for i in 1:n, j in 1:n)
    
    term2 = 0.0
    for i in 1:n
        for j in 1:n
            temp = sum(abs(W[k, i, k, j]) for k in 1:n)
            term2 += abs2(temp)
        end
    end

    term3 = norm(W)
    
    return 0.25 * term1^2 + term2 + 0.25 * term3
end

function find_f_norm_from_w(W, alpha)
    sum1 = norm(W[:, :, :, alpha])
    sum2 = norm(W[:, alpha, :, alpha])
    return 4 * (sum1 - 2* sum2)
end

function simple_norm_operator(W)
    return sum(x -> x^2, W)
end


function calculate_fijkl(alpha, i,j,k,l,eigenvectors,eigenvalues, bonds)
    return (delta(alpha, l)*calculate_Wijkl(i,j,k,alpha,eigenvectors,eigenvalues, bonds)+delta(alpha, k)*calculate_Wijkl(i,j,alpha,l,eigenvectors,eigenvalues, bonds)-delta(alpha, j)*calculate_Wijkl(i,alpha,k,l,eigenvectors,eigenvalues, bonds)-delta(alpha, i)*calculate_Wijkl(alpha,j,k,l,eigenvectors,eigenvalues, bonds))
end

function calculate_fijkl(alpha, i,j,k,l,W)
    return (delta(alpha, l)*W[i,j,k,alpha]+delta(alpha, k)*W[i,j,alpha,l]-delta(alpha, j)*W[i,alpha,k,l]-delta(alpha, i)*W[alpha,j,k,l])
end

function get_f_matrix(alpha, eigenvalues, W)
    F=zeros(Float64,length(eigenvalues),length(eigenvalues),length(eigenvalues),length(eigenvalues))
    for i=1:length(eigenvalues)
        for j=1:length(eigenvalues)
            for k=1:length(eigenvalues)
                for l=1:length(eigenvalues)
                    F[i,j,k,l]=calculate_fijkl(alpha, i,j,k,l,W)
                end
            end
        end
    end
    return F
end

function get_W_separation(eigenvectors, eigenvalues)
    V=build_Vijkl_matrix(eigenvectors)
    energy=build_energy_differences_matrix(eigenvalues)
    W = V./energy
    V_flat=V[:]
    energy_flat=energy[:]
    W_flat=W[:]
    idxs = partialsortperm(W, 1:10, rev=true)
    println("Ws")
    println(W_flat[idxs])
    println("Vs")
    println(V_flat[idxs])
    println("energies")
    println(energy_flat[idxs])
end

function calculate_g_orbital_coefficient(alpha, i, j, k, l,eigenvectors,eigenvalues, bonds, N)
    #println(size(eigenvectors))
    sum=0
    for m in range(1,N)
        for n in range(1,N)
            for o in range(1,N)
                for p in range(1,N)
                    temp=calculate_fijkl(alpha,m,n,o,p,eigenvectors,eigenvalues, bonds)*conj(eigenvectors[:,m][i])*conj(eigenvectors[:,n][j])*eigenvectors[:,o][k]*eigenvectors[:,p][l]
                    sum+=temp
                end
            end
        end
    end
    return sum
end


function calculate_g_norms(alpha, eigenvectors,eigenvalues, bonds, N)  
    sum=0
    count=0
    for i=1:N
        for j=1:N
            for k=1:N
                for l=1:N
                    count+=1
                    println("progress: "*string(count)*"/"*string(N^4))
                    temp=abs2(calculate_g_orbital_coefficient(alpha, i, j, k, l,eigenvectors,eigenvalues, bonds, N))
                    sum+=temp
                end
            end
        end
    end
    return sum
end

function calculate_f_from_w(W, alpha)
    result = zeros(Float64, size(W))

    result1 = zeros(Float64, size(W))

    result2 = zeros(Float64, size(W))

    result3 = zeros(Float64, size(W))
    one=W[:,:,:,alpha]
    two=W[:,:,alpha,:]
    three=W[:,alpha,:,:]
    four=W[alpha,:,:,:]

    result[:,:,:,alpha]=one
    result1[:,:,alpha,:]=two
    result2[:,alpha,:,:]=three
    result3[alpha,:,:,:]=four

    return -result-result1+result2+result3
end

function calculate_g_from_f(F,eigenvectors)
    G= @tensor result[i,j,k,l]:=F[a,b,c,d] * eigenvectors[i,a] * eigenvectors[j, b] * conj.(eigenvectors)[k, c] * conj.(eigenvectors)[l, d]
    return G
end

function find_point(eigenvectors, alpha)
    i=argmax(abs.(eigenvectors[:,alpha]))
    return (i,i,i,i)
end

#q: what is the point?
function get_g_norm_as_function_of_distances(G, eigenvectors, distances, alpha)
    point=find_point(eigenvectors, alpha)
    norms=[]
    for i=1:length(distances)
        distance=distances[i]
        println("distance="*string(distance))
        #indexes=get_indicies_distance(points_matrix, distance, L)
        indexes=generate_coordinates_in_radius(4, distance)
        println("lala")
        indexes=map(tup -> (tup .+ point).%length(eigenvectors[:,1]), indexes)
        #indexes=filter(t -> all(x -> x+1 <= length(eigenvectors[:,1]), t), indexes)
        array=[G[idx.+1...] for idx in indexes]
        push!(norms, norm(array))
    end
    println("??")
    println(length(norms))
    return norms
end

function get_indicies_distance(points, distance, n)
    println("yooo")
    tree = NearestNeighbors.BallTree(points)
    println("hello")
    indices = inrange(tree, points, distance)
    println("hallo")
    return points[indices, :]
end

function get_most_likely(eigenvector)
    return argmax(eigenvector)
end

function to_sum_k_rec(n, k)
    if n == 1
        return ((k,),)
    else
        result = ()
        for x in 0:k-1
            for i in to_sum_k_rec(n - 1, k - x)
                result = ((x,) .+ i,)
            end
        end
        return result
    end
end

#length=how many coordinates
#total_sum=what do you want the coordinates to sum to
function generate_coordinates_in_radius(length, total_sum)
    if length == 1
        return [(total_sum,)]
    else
        result = []
        for value in 0:total_sum
            for permutation in generate_coordinates_in_radius(length - 1, total_sum - value)
                push!(result, (value, permutation...))
            end
        end
        return result
    end
end

function get_rs_of_eigenvectors(eigenvectors)
    rs=[]
    for alpha=1:size(eigenvectors, 2)
        r1=argmax(abs.(eigenvectors[:,alpha]))
        push!(rs, r1)
    end
    return rs
end

function get_distances_matrix(eigenvectors)
    L=size(eigenvectors,2)
    distances=zeros(Float64, (L,L,L,L))

    for i=1:L
        for j=1:L
            for k=1:L
                for l=1:L
                    pos1=argmax(abs.(eigenvectors[:, i]))
                    pos2=argmax(abs.(eigenvectors[:, j]))
                    pos3=argmax(abs.(eigenvectors[:, k]))
                    pos4=argmax(abs.(eigenvectors[:, l]))
                    distance=max(abs(pos1-pos2),abs(pos1-pos3),abs(pos1-pos4),abs(pos2-pos3),abs(pos2-pos4),abs(pos3-pos4))
                    distances[i,j,k,l]=distance
                end
            end
        end
    end
    return distances
end

function transport_operator(eigenvectors, eigenvalues, i_0)
    te=i_0==length(eigenvalues) ? 1 : i_0+1
    result=zeros(Float64, length(eigenvalues), length(eigenvalues))
    Threads.@threads for a=1:length(eigenvalues)
        for b=1:length(eigenvalues)
            if(a==b)
                result[a,b]=0
            else
                result[a,b]=(eigenvectors[i_0, a]*eigenvectors[te, b]+eigenvectors[i_0, b]*eigenvectors[te, a])/(eigenvalues[a]-eigenvalues[b])
            end
        end
    end
    return ifelse.(abs.(result) .< 1e-12, 0.0, result)
end

function interacting_transport_operator(eigenvectors, eigenvalues, i_0)
    O=transport_operator(eigenvectors, eigenvalues, i_0)
    V=Array{Float64, 4}(build_Vijkl_matrix(eigenvectors))
    E=build_energy_differences_matrix(eigenvalues)
    #println(typeof(V))
    #println(typeof(O))
    #println("hello")
    #println(V)
    #println(O)
    result = zeros(Float64, size(V))
    @tensor result[a,b,c,d]:=V[a,b,c,i]*O[i,d]+V[a,b,i,d]*O[i,c]-V[a,i,c,d]*O[i,b]-V[i,b,c,d]*O[i,a]
    #println(result)
    return result./E
end

function interacting_transport_operator_0(eigenvectors, eigenvalues, O)
    #O=transport_operator(eigenvectors, eigenvalues, i_0)
    V=Array{Float64, 4}(build_Vijkl_matrix(eigenvectors))
    E=build_energy_differences_matrix(eigenvalues)
    #println(typeof(V))
    #println(typeof(O))
    #println("hello")
    #println(V)
    #println(O)
    result = zeros(Float64, size(V))
    @tensor result[a,b,c,d]:=V[a,b,c,i]*O[i,d]+V[a,b,i,d]*O[i,c]-V[a,i,c,d]*O[i,b]-V[i,b,c,d]*O[i,a]
    #println(result)
    return result./E
end

#operator is in matrix form
function operator_norm(O)
    temp=transpose(O)*O
    values, __ =eigsolve(temp, 1, :LR; krylovdim=100, ishermitian=true);
    return sqrt(values[1])
end