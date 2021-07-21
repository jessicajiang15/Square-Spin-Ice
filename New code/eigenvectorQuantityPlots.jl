include("hamiltonianDiagonalization.jl")
include("calculations and plotting.jl")



function plotSxVersusH(hList)

end


function findLowest(list)
    min=1;
    for i=2:length(list)
        if(list[i]<list[min])
            min=i;
        end
    end
    return min;
end

function plotSzVersusHTransverse(hList, J, bonds, N)
    szlist=Float64[];
    for h=1:length(hList)
        temp = calculateEigensystemTransverse(N, J, h, bonds,"lanczos", "one", h, 0);
        minInd=findLowest(temp[1]);
        eigenvector = temp[2][minInd];
        #temp[3] is a map mapping index to the list of states
        states= minInd>=temp[3][1][1]&&minInd<=temp[3][1][2] ? temp[3][2] : temp[3][4];
        push!(szlist, calculateSz(eigenvector, states, N));
    end
    return szlist;
end


function generateHListUniform(J, num)
    list::Array{Float64}=Float64[];
    interval=J*10/num;
    i=(1/10)*J;
    while(i<J*10)
        push!(list, i);
        i+=interval;
    end
    return list;
end

function generateHListLog(J, num)
    list::Array{Float64}=Float64[];
    i=0.1*J;
    start=log(2, i);
    count=0;
    while(count<num)
        push!(list,2^start);
        start+=0.01;
        count+=1;
    end
    return list;
end
