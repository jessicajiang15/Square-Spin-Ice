

using LinearAlgebra

function getTi(i, I)
    return mod(I÷(2^(i)), 2);
end


function dir(i, j)
    return i==j ? 1 : 0;
end

#flips the bits at the digits i and j
#credit: https://stackoverflow.com/questions/18247126/how-to-flip-a-bit-at-a-specific-position-in-an-integer-in-any-language/18247246
function flipBits(i, j, n)
    temp=n⊻(1<<i);
    temp=temp⊻(1<<j);
    return temp;
end

function constructMatrix(N, arr)
    #arr contains information about the pairs
    #each arr element is another arr maybe with 3 numbers: J and i and j
    #i is the position of first relevant digit, j second
    tot=2^N;
    H=zeros(Int64, tot, tot)
    for m=1:size(arr, 1)
        #initialize matrix that is totxtot with all 0s
        Htemp=zeros(Int64, tot, tot);
        #loop through all the numbers
        for i=1:tot
            for j=1:tot
                spin1=getTi(arr[m,1], j-1)== 1 ? 1 : -1;
                spin2=getTi(arr[m,2], j-1)== 1 ? 1 : -1;
                if i==j
                    println("i: ",i," j: ",j," spin1spin2: ",spin1*spin2);
                    Htemp[i,j]=spin1*spin2;
                else
                    if(dir(i,flipBits(arr[m,1],arr[m,2],j-1))!=0)
                        temp=1;
                        if(spin1==spin2)
                            temp-=1;
                        else
                            temp+=1;
                        end
                        Htemp[i,j]=temp;
                    end
                end
                end
            end
            H+=arr[m,3]*Htemp;
    end
    return H;

end



function diagonalize(N, arr)
    H=constructMatrix(N, arr);
    println(H);
    eigenvalues=eigvals(H)
    eigenvectors=eigvecs(H)
    temp=Any[];
    push!(temp, eigenvalues);
    push!(temp, eigenvectors);
    return temp;
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
    arr=[0 1 1;1 2 1;2 3 1;3 0 1;0 2 1;1 3 1];
    theResult=diagonalize(N, arr);
    #printEigens(theResult);
end

theRun()
