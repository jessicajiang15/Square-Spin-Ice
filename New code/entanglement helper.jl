include("eigenvectorQuantityPlots.jl")
using TensorCore
using TensorOperations
using ITensors


function getANum(state, listA)
    num=0;
    for i=1:length(listA)
        digit=getTi(listA[i]-1, state);
        num=num|(digit<<(i-1));
    end
    #println(num);
    return num;
end

#e.g. 4, 1 gives 00...0011110
function createMask(digits, pos)
    if(digits<=0)
        #println("wtf");
        return 0;
    end
    num=(2^(digits))-1;
    return num<<pos;
end

#create a number by taking out the bits at positions in listA
function getBNum(state, listA, N)
    #idea: have a mask for each and keep track of last digit that youve finished.
#shift to that digit
#println("state", state);
    num=0;
    lastNum=-1;
    digitsdone=0;
    for i=1:length(listA)
        temp=0;
        numDigits=listA[i]-lastNum-2;
        if(numDigits>0)
            mask=createMask(numDigits, lastNum+1);
            temp=(mask&state)>>(lastNum+1);
            num|=temp<<digitsdone;
            digitsdone+=numDigits;
        end
        lastNum=listA[i]-1;
    end

    numDigits = N - lastNum - 1
    if numDigits > 0
        mask = createMask(numDigits, lastNum + 1)
        temp = (mask & state) >> (lastNum + 1)
        num |= temp << digitsdone
    end
    return num;
end


#create a number by taking out the bits at positions in listA
function getBNumStateNum(state, listA, N)
    #idea: have a mask for each and keep track of last digit that youve finished.
#shift to that digit
#println("state", state);
    num=0;
    lastNum=-1;
    digitsdone=0;
    for i=1:length(listA)
        temp=0;
        numDigits=listA[i]-lastNum-1;
        #println("numdigits",numDigits, "i, ", listA[i]);
        if(numDigits>0)
            mask=createMask(numDigits, lastNum+1);
            #println("mask",mask);
            #println("lastNum",lastNum);
            temp=(mask&state)>>(lastNum+1);
            #println("temp",temp);
            num|=temp<<digitsdone;
            #println("num",temp);
            digitsdone+=numDigits;
            #println("digits done ",digitsdone);
        end
        lastNum=listA[i];
    end

    numDigits = N - lastNum - 1
    if numDigits > 0
        mask = createMask(numDigits, lastNum + 1)
        temp = (mask & state) >> (lastNum + 1)
        num |= temp << digitsdone
    end

    return num
end




function constructCoefficientMatrix(eigenvector,states,listA, N)
    #listA is a list of the positions of interest
    H::Matrix{Complex}=zeros(Complex, 2^length(listA), 2^(N-length(listA)));
    #zeros(Float64, length(list),length(list));
    for i=1:length(states)
        #println("listA", listA)
        #1001
        #println("lallalalala")
        getA=getANum(states[i], listA);
        getB=getBNum(states[i], listA, N);

        #println("A: "*string(getA))
        #println("B: "*string(getB))

        H[getA+1,getB+1]=eigenvector[i];
    end
    #println(H);
    return H;
end

function getEntanglementEntropy(eigenvector, states, listA, N)
    H=constructCoefficientMatrix(eigenvector, states, listA, N);
    decomp=svd(H);
    tot=0;
    for i=1:length(decomp.S)
        temp=decomp.S[i];
        #println(temp);
        tot+= temp == 0 ? 0 : -(temp^2)*log(temp^2);
    end
    return tot;
end

#truncate up to n
function getApproxEnergy(eigenvector, states, listA, N, k)
    H=constructCoefficientMatrix(eigenvector, states, listA, N);
    decomp=svd(H);


    ss=decomp.S[1:k]
    uu=decomp.U[:,1:k]
    vv=decomp.V[:,1:k]

    test=(uu)*Diagonal(ss)*Transpose(vv)

    vector=deconstructCoefficientMatrix(test,states,listA, N)
#=
    vector=zeros(size(uu, 1) * size(vv, 1))
    for i=1:k
        vector=vector+ss[i]*kron(uu[:,i],vv[:,i])
    end
=#
    error=0
    for i=k+1:length(decomp.S)
        error+=decomp.S[i]^2
    end

    return (error,vector/norm(vector))
end

function getApproximateEnergy(eigenvector, states, listA, N, k)
    H=getTruncatedHamiltonian(eigenvector, states, listA, N, k)
    
    return uu*Digaonal(ss[1:k])*vv
end

function matrix_product_state(eigenvector, N, k)
    listA=[1]
    eigvec=copy(eigenvector)
    states=[]

    for i=0:2^N-1
        push!(states, i)
    end

    H=constructCoefficientMatrix(eigenvector, states, listA, N);
    println(size(H))
    A_list=[]
    ks=[]
    for i=1:N-1
        decomp=svd(H);
        kk=min(k, size(H, 1),size(H, 2))
        #println(size(H, 1))
        push!(ks, kk)

        A_1=decomp.U[:,1:kk]

        #println(decomp.U[:,1:kk])


        A_2=Diagonal(decomp.S[1:kk])*Transpose(decomp.V[:,1:kk])

        #A_2=reshape(A_2, i+1, k, N)
        A_2=Base.reshape(A_2, (2*kk, 2^(N-i-1)))
        H=A_2

        if i!=N-1
            if(i!=1)
                A_1=reshape_3d(A_1,ks[end-1],kk)
                #reshape_3d(Transpose(A_1), ks[end-1],kk)
            end
            push!(A_list, A_1)
        else
            A_1=reshape_3d(A_1,ks[end-1],kk)
            #reshape_3d((A_1), ks[end-1],kk)
            push!(A_list, A_1)
            push!(A_list, Transpose(Diagonal(decomp.S[1:kk])*Transpose(decomp.V[:,1:kk])))
        end
    end
    
    return A_list
end

# l is which iteration this is– 0,1,2,3...
# matrix product state helper
function reshape_1(H, l, k, N)

    H1::Matrix{Float64}=zeros(Float64, k*2, 2^(N-l));
    listA=[1]
    states=[]
    for j=1:k
        println(N-l+1)
        for i=0:2^(N-l+1)-1
            getA=getANum(i, listA);
            getB=getBNum(i, listA, N);
            H1[2*(k-1)+getA+1, getB+1]=H[j, i+1]
        end
    end
    return H1
end

function reshape_3d(H, k_prev,k)
    H1=zeros(Float64, 2, k_prev, k);

    for i=0:k_prev-1
        for j=0:k-1
            for l=1:2
                H1[l,i+1,j+1]=H[(l-1)*k_prev+i+1,j+1]
            end
        end
    end
    return H1
end

function reshape_3d_new(H, k_prev,k)
    H1=zeros(Float64, 2, k_prev, k);

    for i=0:k_prev-1
        for j=0:k-1
            for l=1:2
                #H1[l,i+1,j+1]=H[(j)*2+l,i+1]
                H1[l,i+1,j+1]=H[(l-1)*k+j+1,i+1]
            end
        end
    end
    return H1
end

#assume real indicies are first one, contract them first
function calculate_mps_overlap(a_list_1, a_list_2)
    contracted_list=[]
    for i=1:(length(a_list_1))
        if(i!=1&&i!=length(a_list_1))
        @tensor new_result[b,c,d,e] := a_list_1[i][a, b, c] * a_list_2[i][a, d, e]
        push!(contracted_list,new_result)
        else
        @tensor new_result[b,c] := a_list_1[i][a, b] * a_list_2[i][a, c]
        push!(contracted_list,new_result)
        end
    end
    temp=[contracted_list[1]];

    for i=2:length(contracted_list)-1
        @tensor new_result[c, d]:=temp[end][a, b]*contracted_list[i][a, c, b, d]
        push!(temp, new_result)
    end
   # @tensor new_result := sum(temp[end][a, b]*contracted_list[end][a,b])
return tr(temp[end] * contracted_list[end])
end

function construct_H_tensor(H, L)
    num_dims=2*L
    dims = fill(2, num_dims)
    tensor = zeros(dims...)

    for i=0:size(H,0)-1
        i_index=[]
        for k=1:L
            push!(indexes, getTi(k,i)+1)
        end
        for j=0:size(H,1)-1
            ii_index=copy(i_index)
            j_index=[]
            for k=1:L
                push!(j_index, getTi(k,j)+1)
            end
            append!(ii_index, j_index)
            tensor[Tuple(ii_index)]=H[i+1, j+1]
        end
    end
    return tensor
end

function apply_sx_then_contract(ii, a_list_1, a_list_2, h)
    sx=[0 1;1 0]
    contracted_list=[]
    for i=1:(length(a_list_1))
        if(i!=1&&i!=length(a_list_1))
        if(i==ii)
        @tensor new_result[b,c,d,e] := a_list_1[i][a, b, c] *sx[f,a]* a_list_2[i][f, d, e]
        else
        @tensor new_result[b,c,d,e] := a_list_1[i][a, b, c]* a_list_2[i][a, d, e]
        end
        push!(contracted_list,new_result)
        else
        if(i==ii)
            @tensor new_result[b,c] := a_list_1[i][a, b] * sx[a,f]*a_list_2[i][f, c]
        else
            @tensor new_result[b,c] := a_list_1[i][a, b] * a_list_2[i][a, c]
        end
        push!(contracted_list,new_result)
        end
    end
    temp=[contracted_list[1]];

    for i=2:length(contracted_list)-1
        @tensor new_result[c, d]:=temp[end][a, b]*contracted_list[i][a, c, b, d]
        push!(temp, new_result)
    end
   # @tensor new_result := sum(temp[end][a, b]*contracted_list[end][a,b])
return -h*tr(temp[end] * contracted_list[end])
end



function apply_sz_then_contract(ii, a_list_1, a_list_2, h)
    sz=[1 0;0 -1]
    contracted_list=[]
    for i=1:(length(a_list_1))
        if(i!=1&&i!=length(a_list_1))
        if(i==ii)
        @tensor new_result[b,c,d,e] := a_list_1[i][a, b, c] *sz[f,a]* a_list_2[i][f, d, e]
        else
        @tensor new_result[b,c,d,e] := a_list_1[i][a, b, c]* a_list_2[i][a, d, e]
        end
        push!(contracted_list,new_result)
        else
        if(i==ii)
            @tensor new_result[b,c] := a_list_1[i][a, b] * sz[a,f]*a_list_2[i][f, c]
        else
            @tensor new_result[b,c] := a_list_1[i][a, b] * a_list_2[i][a, c]
        end
        push!(contracted_list,new_result)
        end
    end
    temp=[contracted_list[1]];

    for i=2:length(contracted_list)-1
        @tensor new_result[c, d]:=temp[end][a, b]*contracted_list[i][a, c, b, d]
        push!(temp, new_result)
    end
   # @tensor new_result := sum(temp[end][a, b]*contracted_list[end][a,b])
return -h*tr(temp[end] * contracted_list[end])
end


function apply_sz_then_contract(ii, jj, a_list_1, a_list_2, J)
    sz=[1 0;0 -1]
    contracted_list=[]
    for i=1:(length(a_list_1))
        if(i!=1&&i!=length(a_list_1))
        if(i==ii||i==jj)
            @tensor new_result[b,c,d,e] := a_list_1[i][a, b, c] * sz[a, f] * a_list_2[i][f, d, e]
        else
            @tensor new_result[b,c,d,e] := a_list_1[i][a, b, c] * a_list_2[i][a, d, e]
        end
        push!(contracted_list,new_result)
        else
            if(i==ii||i==jj)
                @tensor new_result[b,c] := a_list_1[i][a, b] * sz[a, f]*a_list_2[i][f, c]
            else
                @tensor new_result[b,c] := a_list_1[i][a, b] * a_list_2[i][a, c]
            end
        push!(contracted_list,new_result)
        end
    end
    temp=[contracted_list[1]];

    for i=2:length(contracted_list)-1
        @tensor new_result[c, d]:=temp[end][a, b]*contracted_list[i][a, c, b, d]
        push!(temp, new_result)
    end
   # @tensor new_result := sum(temp[end][a, b]*contracted_list[end][a,b])
return J*tr(temp[end] * contracted_list[end])
end

function get_energy_mps(a_list_1, a_list_2, bonds, J, h)
    sum=0
    for i=1:length(a_list_1)
       sx=apply_sx_then_contract(i, a_list_1, a_list_2, h)
       sum+=sx
    end
    temp=0
    for b in bonds
        bond1::Int=b.site1.num;
        bond2::Int=b.site2.num;
        sz=apply_sz_then_contract(bond1, bond2, a_list_1, a_list_2, J)
        temp+=sz
        sum+=sz
    end
    return sum
end

#14 sites
function deconstruct_mps(a_list, N)
    result = a_list[1]
    indices = ["s1", "v1"]
    new_indices = ["s1"]
@tensor result[s1,s2,v2]:=a_list[1][s1,v1]*a_list[2][s2,v1,v2]
@tensor result[s1,s2,s3,v3]:=result[s1,s2,v2]*a_list[3][s3,v2,v3]
@tensor result[s1,s2,s3,s4,v4]:=result[s1,s2,s3,v3]*a_list[4][s4,v3,v4]
@tensor result[s1,s2,s3,s4,s5,v5]:=result[s1,s2,s3,s4,v4]*a_list[5][s5,v4,v5]
@tensor result[s1,s2,s3,s4,s5,s6,v6]:=result[s1,s2,s3,s4,s5,v5]*a_list[6][s6,v5,v6]
@tensor result[s1,s2,s3,s4,s5,s6,s7,v7]:=result[s1,s2,s3,s4,s5,s6,v6]*a_list[7][s7,v6,v7]
@tensor result[s1,s2,s3,s4,s5,s6,s7,s8,v8]:=result[s1,s2,s3,s4,s5,s6,s7,v7]*a_list[8][s8,v7,v8]
@tensor result[s1,s2,s3,s4,s5,s6,s7,s8,s9,v9]:=result[s1,s2,s3,s4,s5,s6,s7,s8,v8]*a_list[9][s9,v8,v9]
@tensor result[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,v10]:=result[s1,s2,s3,s4,s5,s6,s7,s8,s9,v9]*a_list[10][s10,v9,v10]
@tensor result[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,v11]:=result[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,v10]*a_list[11][s11,v10,v11]
@tensor result[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,v12]:=result[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,v11]*a_list[12][s12,v11,v12]
@tensor result[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,v13]:=result[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,v12]*a_list[13][s13,v12,v13]
@tensor result[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14]:=result[s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,v13]*a_list[14][s14,v13]

println(size(result))
    # Assuming N is defined somewhere in your scope
    return reshape(result,2^N)
end


function deconstruct_mps_2(a_list_1)
    a_list_temp=[]
    s_Strings=[]

    s="s1"
    push!(s_Strings, s)
    i_String="i1"
    temp_size=size(a_list_1[1])

    s1=Index(temp_size[1],s)
    ii=Index(temp_size[2],i_String)
    push!(a_list_temp,ITensor(a_list_1[1], s1, ii))
    
    for i=2:length(a_list_1)-1
        s="s"*string(i)
        push!(s_Strings, s)

        i_String="i"*string(i)
        i_String_1="i"*string(i-1)
            temp_size=size(a_list_1[i])
            s1=Index(temp_size[1],s)
            ii=Index(temp_size[2],i_String_1)
            ii1=Index(temp_size[3],i_String)
            push!(a_list_temp,ITensor(a_list_1[i], s1, ii, ii1))
    end

    s="s"*string(length(a_list_1))
    push!(s_Strings, s)

    i_String="i"*string(length(a_list_1)-1)
    temp_size=size(a_list_1[end])

    s1=Index(temp_size[1],s)
    ii=Index(temp_size[2],i_String)
    push!(a_list_temp,ITensor(a_list_1[end], s1, ii))

    result=a_list_temp[1]
    for i=2:length(a_list_temp)
        result=result*a_list_temp[i]
    end

    #println(result)
    result_1=Array(result)

    return destroy(result_1)
end

function destroy(tensor)
    N=ndims(tensor)

    H1::Vector{Float64}=zeros(Float64, 2^length(tensor)-1);
    for i=0:2^length(tensor)-1
        position=[]
        for j=1:N
            index=getTi(j, i)
            push!(position, index+1)
        end
        H1[i+1]=tensor[Tuple(position)]
    end
    return H1
end

function deconstructCoefficientMatrix(H,states,listA, N)
    #listA is a list of the positions of interest
    eigenvector=zeros(Float64, 2^N);
    #zeros(Float64, length(list),length(list));
    for i=1:length(states)
        #println("listA", listA)
        #1001
        getA=getANum(states[i], listA);
        getB=getBNum(states[i], listA, N);
        #println("A", getA)
        #println("B", getB)

        eigenvector[i]=H[getA+1,getB+1];
    end
    #println(H);
    return eigenvector;
end


function calculate_correlation_mps(a_list,r,N)
    sz_1=apply_sz_then_contract(1, r, a_list, a_list, 1)
    return sz_1
end

function change_to_canoncial_form(a_list, l,k)
    new_a_list=[]
    W=[]
    values=[]
    @tensor W[s1, s2, v2]:=a_list[1][s1,v1]*a_list[2][s2,v1,v2]
    W=reshape(W, (2, 2*size(W, 3)))
    temp=svd(W)
    push!(new_a_list, temp.U[:,1:size(a_list[1],2)])
    a_list[2]=reshape_3d_new(copy(Transpose(Diagonal(temp.S)*Transpose(temp.V))), size(a_list[2],2), size(a_list[2],3))
    #a_list[2]=reshape_3d(Transpose(Diagonal(temp.S)*Transpose(temp.V)), size(a_list[2],2), size(a_list[2],3))
    for i in range(2, length(a_list)-1)
        if(i==length(a_list)-1)
            @tensor W[s1,v1,s2]:=a_list[i][s1,v1, v2]*a_list[i+1][s2,v2]            
            W=reshape(W, 2*size(W, 2), 2)
            temp=svd(W)
            push!(new_a_list, reshape_3d(temp.U, size(a_list[end-1], 2), size(a_list[end-1], 3)))
            push!(new_a_list, Transpose(Diagonal(temp.S[1:size(a_list[end],2)])*Transpose(temp.V[:,1:size(a_list[end],2)])))
        else
        @tensor W[s1, v1,s2,v3]:=a_list[i][s1,v1,v2]*a_list[i+1][s2,v2,v3]
        W=reshape(W, (2*size(a_list[i], 2), 2*size(a_list[i+1], 3)))

        temp=svd(W)
        push!(new_a_list, reshape_3d(temp.U, size(a_list[i], 2), size(a_list[i], 3)))
        #a_list[i+1]=reshape_3d(Transpose(Diagonal(temp.S)*Transpose(temp.V)), size(a_list[i+1],2), size(a_list[i+1],3))
        println("i "*string(i))
        println(size(a_list[i+1],2))
        a_list[i+1]=reshape_3d_new(Transpose(Diagonal(temp.S[1:(size(a_list[i+1],2))])*Transpose(temp.V[:,1:(size(a_list[i+1],2))])), size(a_list[i+1],2), size(a_list[i+1],3))
        #a_list[i+1]=reshape(copy(Transpose(Diagonal(temp.S[1:min(size(a_list[i+1],2))])*Transpose(temp.V[:,1:min(size(a_list[i+1],2))]))), (2, size(a_list[i+1],2), size(a_list[i+1],3)))
        end
    end

    println(size(new_a_list[end-1],3))
    println(size(new_a_list[end],2))

    @tensor W[s1, v1, s2]:=new_a_list[end-1][s1,v1, v2]*new_a_list[end][s2,v2]
    W=reshape(W, (2*size(W, 2), 2))
    temp=svd(W)
    new_a_list[end]=((temp.V[:,1:size(a_list[end],2)]))
    push!(values, temp.S[1:2])
    new_a_list[end-1]=reshape_3d(temp.U*Diagonal(temp.S), size(a_list[end-1], 2), size(a_list[end-1], 3))

    for i=1:length(new_a_list)
        println(size(new_a_list[i]))
    end

    for i in reverse(l:length(a_list)-1)
        if (i!==2)
            @tensor W[s1, v1,  s2, v3]:=new_a_list[i-1][s1,v1, v2]*new_a_list[i][s2,v2, v3]
            W=reshape(W, (2*size(W, 2), 2*size(W, 4)))
            kk=k
            #size(W, 1),size(W, 2)

            temp=svd(W)
            new_a_list[i]=reshape_3d_new(copy(Transpose(Transpose(temp.V[:,1:min(kk,size(a_list[i], 2))]))), min(kk,size(a_list[i], 2)), min(kk,size(a_list[i], 3)))
            #new_a_list[i]=reshape_3d(Transpose(Transpose(temp.V[:,1:kk])), size(a_list[i], 2), size(a_list[i], 3))
            push!(values, temp.S[1:2])
            index=min(kk,size(a_list[i], 2),size(a_list[i-1], 3))
            new_a_list[i-1]=reshape_3d(copy(temp.U[:,1:index]*Diagonal(temp.S[1:index])), size(a_list[i-1], 2), min(kk,size(a_list[i-1], 3)))
            #new_a_list[i-1]=reshape_3d(temp.U[:,1:kk]*Diagonal(temp.S[1:kk]), size(a_list[i-1], 2), size(a_list[i-1], 3))
        else
            @tensor W[s1, s2, v2]:=new_a_list[1][s1,v1]*new_a_list[2][s2,v1, v2]
            W=reshape(W, (2, 2*size(a_list[i], 3)))
            temp=svd(W)
            kk=min(k, size(W, 1),size(W, 2))
            new_a_list[2]=reshape_3d_new(copy(Transpose(Transpose(temp.V[:,1:min(kk,size(a_list[i], 2))]))), min(kk,size(a_list[i], 2)), size(a_list[i], 3))
            #new_a_list[2]=reshape_3d(Transpose(Transpose(temp.V[:,1:min(kk,size(a_list[i], 2))])), min(kk,size(a_list[i], 2)), size(a_list[i], 3))
            #new_a_list[1]=temp.U[:,1:min(kk,size(a_list[i], 2))]*Diagonal(temp.S[1:min(kk,size(a_list[i], 2))])
        end
    end
    return new_a_list
end



function change_to_canoncial_form(a_list, l)
    new_a_list=[]
    W=[]
    values=[]
    @tensor W[s1, s2, v2]:=a_list[1][s1,v1]*a_list[2][s2,v1,v2]
    W=reshape(W, (2, 2*size(W, 3)))
    temp=svd(W)
    push!(new_a_list, temp.U)
    a_list[2]=reshape_3d(Transpose(Diagonal(temp.S)*Transpose(temp.V)), 2, length(temp.S))
    for i in range(2, length(a_list)-1)
        if(i==length(a_list)-1)
            @tensor W[s1,v1,s2]:=a_list[i][s1,v1, v2]*a_list[i+1][s2,v2]
            W=reshape(W, 2*size(W, 2), 2)
            temp=svd(W)
            push!(new_a_list, reshape_3d(temp.U, size(a_list[end-1], 2), size(a_list[end-1], 3)))
            push!(new_a_list, Transpose(Diagonal(temp.S)*Transpose(temp.V)))
        else
        @tensor W[s1, v1,s2,v3]:=a_list[i][s1,v1,v2]*a_list[i+1][s2,v2,v3]
        W=reshape(W, (2*size(a_list[i], 2), 2*size(a_list[i+1], 3)))
        temp=svd(W)
        push!(new_a_list, reshape_3d(temp.U, size(a_list[i], 2), size(a_list[i], 3)))
        a_list[i+1]=reshape_3d(Transpose(Diagonal(temp.S)*Transpose(temp.V)), size(a_list[i+1],2), size(a_list[i+1],3))
        end
    end


    @tensor W[s1, v1, s2]:=new_a_list[end-1][s1,v1, v2]*new_a_list[end][s2,v2]
    W=reshape(W, (2*size(W, 3), 2))
    temp=svd(W)
    new_a_list[end]=(Transpose(temp.V))
    new_a_list[end-1]=reshape_3d(temp.U*Diagonal(temp.S), size(a_list[end-1], 2), size(a_list[end-1], 3))

    for i in reverse(l:length(a_list)-1)
        if (i!==2)

            @tensor W[s1,v1,  s2, v3]:=new_a_list[i-1][s1,v1, v2]*new_a_list[i][s2,v2, v3]
            W=reshape(W, (2*size(a_list[i-1], 2), 2*size(a_list[i], 3)))
            temp=svd(W)
            new_a_list[i]=reshape_3d(Transpose(Transpose(temp.V)), size(a_list[i], 2), size(a_list[i], 3))
            push!(values, temp.S[1:2])
            new_a_list[i-1]=reshape_3d(temp.U*Diagonal(temp.S), size(a_list[i-1], 2), size(a_list[i-1], 3))
        else
            @tensor W[s1, s2, v2]:=new_a_list[1][s1,v1]*new_a_list[2][s2,v1, v2]
            W=reshape(W, (2, 2*size(a_list[i], 3)))
            temp=svd(W)
            new_a_list[2]=reshape_3d(Transpose(Transpose(temp.V)), size(a_list[i], 2), size(a_list[i], 3))
            push!(values, temp.S[1:2])
            new_a_list[1]=temp.U*Diagonal(temp.S)
        end
    end
    return (new_a_list, values)
end

function apply_trotter_gate(a_list, i, j, deltat)
    sz=[1 0 0 0 ; 0 -1 0 0 ; 0 0 -1 0 ; 0 0 0 1]
    matrix=exp(-sz*deltat)
    matrix=reshape(matrix, (2,2,2,2))

    a_1=a_list[i]
    a_2=a_list[j]
    #@tensor new_a[s1,k1,s2,k3]:=a_1[s1, k1, k2]*a_2[s2, k2, k3]
    #@tensor W[s1,v1,  s2, v3]:=new_a_list[i-1][s1,v1, v2]*new_a_list[i][s2,v2, v3]
    @tensor new_a[s1, k1, s2, k3] := a_1[s1, k1, k2] * a_2[s2, k2, k3]
    @tensor new_a[d1,k1,d2,k3]:=matrix[s1,s2,d1,d2]*new_a[s1,k1,s2,k3]
    return new_a
end

function apply_trotter_gate_begin(a_list, deltat)
    sz=[1 0 0 0 ; 0 -1 0 0 ; 0 0 -1 0 ; 0 0 0 1]
    matrix=exp(-sz*deltat)
    matrix=reshape(matrix, (2,2,2,2))
    a_1=a_list[1]
    a_2=a_list[2]
    @tensor new_a[s1,s2,k2]:=a_1[s1, k1]*a_2[s2, k1, k2]
    @tensor new_a[d1,d2,k2]:=matrix[s1,s2,d1,d2]*new_a[s1,s2,k2]
    return new_a
end

function apply_trotter_gate_end(a_list, deltat)
    sz=[1 0 0 0 ; 0 -1 0 0 ; 0 0 -1 0 ; 0 0 0 1]
    matrix=exp(-sz*deltat)
    matrix=reshape(matrix, (2,2,2,2))
    a_1=a_list[end-1]
    a_2=a_list[end]
    @tensor new_a[s1,k1,s2]:=a_1[s1, k1,k2]*a_2[s2, k2]
    @tensor new_a[d1,k1,d2]:=matrix[s1,s2,d1,d2]*new_a[s1,k1,s2]
    return new_a
end

function apply_trotter_gate_field(a_list, i, deltat, hx, hz)
    szsx=[-hz -hx ; -hx hz]
    matrix=exp(-szsx*deltat)
    a_1=a_list[i]
    @tensor new_a[d1,k1,k2]:=matrix[d1,s1]*a_1[s1, k1, k2]
    return new_a
end

function apply_trotter_gate_field_begin(a_list,deltat, hx, hz)
    szsx=[-hz -hx ; -hx hz]
    matrix=exp(-szsx*deltat)
    a_1=a_list[1]
    @tensor new_a[d1,k1]:=matrix[d1,s1]*a_1[s1, k1]
    return new_a
end

function apply_trotter_gate_field_end(a_list,deltat, hx, hz)
    szsx=[-hz -hx ; -hx hz]
    matrix=exp(-szsx*deltat)
    a_1=a_list[end]
    @tensor new_a[d1,k1]:=matrix[d1,s1]*a_1[s1, k1]
    return new_a
end


function get_even_bonds(bonds)
    bonds_new=[]
    for bond in bonds
        if(bond.site1.num % 2==0)
            push!(bonds_new, bond)
        end
    end
    return bonds_new
end

function get_odd_bonds(bonds)
    bonds_new=[]
    for bond in bonds
        if(bond.site1.num % 2!=0)
            push!(bonds_new, bond)
        end
    end
    return bonds_new
end

function apply_trotter_gates_even(a_list, deltat, bonds)
    even_bonds=get_even_bonds(bonds)
    new_as=[]
    for bond in even_bonds
        i=bond.site1.num
        j=bond.site2.num
        if(i!=1&&i!=length(a_list)-1)
            temp=apply_trotter_gate(a_list, i, j, deltat)
        elseif(i==1)
            temp=apply_trotter_gate_begin(a_list, deltat)
        else
            temp=apply_trotter_gate_end(a_list, deltat)
        end
        push!(new_as, temp)
    end
    return new_as
end

function apply_trotter_gates_field(a_list, deltat, bonds, hx, hz)
    new_as=[]
    for i=1:length(a_list)
        if(i==1)
            temp=apply_trotter_gate_field_begin(a_list, deltat, hx, hz)
        elseif(i==length(a_list))
            temp=apply_trotter_gate_field_end(a_list, deltat, hx, hz)
        else
            temp=apply_trotter_gate_field(a_list, i, deltat, hx, hz)
        end
        push!(new_as,temp)
    end
    return new_as
end

function apply_trotter_gates_odd(a_list, deltat, bonds)
    even_bonds=get_odd_bonds(bonds)
    new_as=[]
    for bond in even_bonds
        i=bond.site1.num
        j=bond.site2.num
        if(i==1)
            temp=apply_trotter_gate_begin(a_list, deltat)
        elseif(i==length(a_list)-1)
            temp=apply_trotter_gate_end(a_list, deltat)
        else
            temp=apply_trotter_gate(a_list, i, j, deltat)
        end
        push!(new_as, temp)
    end
    return new_as
end

#each is a rank 3 or 4 tensor rather than rank 2 or 3 tensor
# assume even

function end_tensor_svd(a, a_list)
    W=reshape(a, (2*size(a_list[end-1], 2), 2))
    temp=svd(W)
    left=reshape_3d(temp.U, size(a_list[end-1], 2), 2*size(a_list[end], 3))
    right=Transpose(Diagonal(temp.S)*Transpose(temp.V))
    return [left, right]
end

function beginning_tensor_svd(a, a_list)
    W=reshape(a, (2, 2*size(a_list[2], 3)))
    temp=svd(W)
    left=temp.U
    right=reshape_3d_new(Transpose(Diagonal(temp.S)*Transpose(temp.V)), size(temp.U,2),size(a_list[2],3))
    #right=reshape_3d(Diagonal(temp.S)*Transpose(temp.V),2*size(a_list[2],2),size(a_list[2],3))
    return [left, right]
end

function middle_tensor_svd(a, a_list, i)
    W=reshape(a, (2*size(a, 2), 2*size(a, 4)))
    temp=svd(W)
    left=reshape_3d(temp.U, size(a,2),size(temp.U,2))
    #reshape_3d(temp.U, size(a_list[i],2),2*size(a_list[i],3))
    right=reshape_3d_new(Transpose(Diagonal(temp.S)*Transpose(temp.V)), size(temp.U,2),size(a_list[i+1],3))
    #right=reshape(Diagonal(temp.S)*Transpose(temp.V),(2,2*size(a_list[i+1],2),size(a_list[i+1])))
    #right=reshape_3d(Diagonal(temp.S)*Transpose(temp.V),2*size(a_list[i+1],2),size(a_list[i+1],3))
    return [left, right]
end

function restore_mps_form(a_list_bad, a_list, bonds,even)
    new_a_list=[]
    if(even)
        push!(new_a_list, a_list[1])
        count=2
        for i=1:length(a_list_bad)
            append!(new_a_list, middle_tensor_svd(a_list_bad[i], a_list, count))
            count+=2
        end
        push!(new_a_list, a_list[end])
    else
        append!(new_a_list, beginning_tensor_svd(a_list_bad[1], a_list))
        count=3
        for i=2:length(a_list_bad)-1
            append!(new_a_list, middle_tensor_svd(a_list_bad[i], a_list, count))
            count+=2
        end
        append!(new_a_list, end_tensor_svd(a_list_bad[end], a_list))
    end
    return new_a_list
end

function time_evolve_step(a_list,deltat, bonds, k, hx, hz)
    a_list_temp=apply_trotter_gates_field(a_list,deltat, bonds, hx, hz)
    a_list_even=apply_trotter_gates_even(a_list_temp, deltat, bonds)
    a_list_temp=restore_mps_form(a_list_even, a_list, bonds,true)

    #a_list_odd=apply_trotter_gates_odd(a_list_temp, deltat, bonds)
    #a_list_temp=restore_mps_form(a_list_odd, a_list_temp, bonds,false)
    return (change_to_canoncial_form(a_list_temp, 2, k))
end

function get_energy_mps_2(a_list_1, a_list_2, bonds, J, h, h2)
    sum=0
    for i=1:length(a_list_1)
       sx=apply_sx_then_contract(i, a_list_1, a_list_2, h)
       sum+=sx
    end

    for i=1:length(a_list_1)
        sz=apply_sz_then_contract(i, a_list_1, a_list_2, h2)
        sum+=sz
     end

    temp=0
    for b in bonds
        bond1::Int=b.site1.num;
        bond2::Int=b.site2.num;
        sz=apply_sz_then_contract(bond1, bond2, a_list_1, a_list_2, J)
        temp+=sz
        sum+=sz
    end
    return sum
end