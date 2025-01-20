using LinearAlgebra
using SparseArrays
using Arpack
using KrylovKit
using Plots

#flips the bits at the digits i and j
#credit: https://stackoverflow.com/questions/18247126/how-to-flip-a-bit-at-a-specific-position-in-an-integer-in-any-language/18247246
function flipBits(i, j, n)
    temp=n⊻(1<<i);
    temp=temp⊻(1<<j);
    return temp;
end

#gets the ith binary digit (Starting from 0) of the decimal number I
function getTi(i, I)
    return mod(I÷(2^(i)), 2);
end

#constructs the hamiltonian. code inspired by "Computational Quantum Spin Systems" by Sandvik and thanks to Yutan for his suggestions and guidance!
function constructMatrix(N, J)
    tot::Int=2^N;
    H::SparseMatrixCSC{Float64}=spzeros(Float64, tot, tot);
    #loop over all the states
        for i=0:tot-1
            #bonds are numbered from 0 to 1
            for j=0:N-1
                #loop over all bonds, avoid double counting by starting from j+1
                #if we do this
                for temp=j+1:N-1
                    if getTi(j, i)==getTi(temp, i)
                        #match=product of spins=1 for z terms, which are the diagonal entries
                        H[i+1, i+1]+=J/4;
                    else
                        #if they don't match the product of the spins is -1 for the sigma z term
                        H[i+1, i+1]-=J/4;
                        #there will be a nonzero entry at the column numbered with the flipped bit version of i
                        #due to orthonormality
                        b=flipBits(temp, j, i);
                        H[i+1, b+1]=J/2;
                    end
                end
            end
        end
        return H;
    end


function diagonalize(N, J)
    H::SparseMatrixCSC{Float64}=constructMatrix(N, J)
    println(H);
    #using eigs from Arpack
    eig=eigs(H, nev=2^N);
    return eig;
end

function printEigens(arr)
    println("Eigenvalues:");
    for i=1:length(arr[1])
        println(arr[1][i])
    end
    println("Eigenvectors:");
    for i=1:length(arr[2])
        println(arr[2][i])
    end


end

function theRun()
    N=4;
    J=1;
    theResult=diagonalize(N, J);

    println(theResult);
    #printEigens(theResult);
end

theRun()
