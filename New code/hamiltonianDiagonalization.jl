include("hamiltonianGeneration.jl")
using LinearAlgebra
using KrylovKit
using Arpack

function singleOutNUpSpins(N, range)
    #temp2 stores the map and the list
    temp2::Array{Any}=Any[];
    #temp stores the statess
    temp::Array{Int}=Int[];
    #stores the states and their corresponding number
    temp1::Dict{Int,Int}=Dict{Int, Int}();
    in::Int=0;
    for i=0:range-1
        count::Int=0;
        n::Int=i;
        while n>0
            count+=1;
            if count>N
                break;
            end
            n=n & (n-1);
        end
        if(count==N)
            in+=1;
            temp1[i]=in;
            push!(temp, i)

    end
end
println("THE STATES",length(temp));
push!(temp2, temp);
push!(temp2, temp1);
    return temp2;

end

function calculateEigensystemHeisenberg(N, J, bonds)
        eigensystem::Array{Any}=Any[];

    eigenvalues::Array{Any}=Any[];
    eigenvectors::Array{Any}=Any[];
    for i=0:N*N
        println("next:");
        spinUps::Array{Any}=singleOutNUpSpins(i, 2^(N*N));
        Htemp::SparseMatrixCSC{Float64}=constructHeisenbergHamiltonian(spinUps, bonds, N, J);
        #    Htemp::Matrix{Float64}=constructHeisenbergHamiltonian(spinUps, bonds, N, J);
        #println(Htemp);
        #TODO:check this, how do you incorporate eigenvectors
        if i==0||i==N*N
            append!(eigenvalues, Htemp[1, 1]);
            local eigenvector::Array{Vector{Float64}}=Array{Vector{Float64}}(undef, 1);
            eigenvector[1,1]=Float64[];
            push!(eigenvector[1,1],1);
            append!(eigenvectors, eigenvector)
            continue;
        end
        eigtemp=eigs(Symmetric(Htemp), nev=length(spinUps[1]));
        #eigtemp=eigen(Htemp);
        #println(eigtemp);
        println(length(eigtemp[1]));

        append!(eigenvalues, eigtemp[1]);
        append!(eigenvectors, eigtemp[2]);

        #append!(eigenvalues, eigtemp.values);
        #append!(eigenvectors, eigtemp.vectors);
    end

push!(eigensystem, eigenvalues);
push!(eigensystem, eigenvectors);

return eigensystem;

end



function calculateEigensystemTransverse(N, J, h, bonds)
    println("timing matrix generation");
    @time begin
        eigensystem::Array{Any}=Any[];
    eigenvalues::Array{Any}=Any[];
    eigenvectors::Array{Any}=Any[];
    evenSpins=singleOutEvenOddSpins(true, 2^(N*N), N);
        oddSpins=singleOutEvenOddSpins(false, 2^(N*N), N);
        println("STARTING EVEN");
        HtempEven=constructTransverseHamiltonian(evenSpins, bonds, N, J,h);
        #HtempEven::Matrix{Float64}=constructHamiltonianTransverse(evenSpins, bonds, N, J,h);
        #println(HtempEven);
        #println(HtempOdd);

        eigtemp=eigs(Symmetric(HtempEven), nev=2^(N*N-1));
        HtempEven=nothing;
        evenSpins=nothing;
        GC.gc();
        println("STARTING ODD");
        HtempOdd=constructTransverseHamiltonian(oddSpins, bonds, N, J, h);
                #HtempOdd::Matrix{Float64}=constructHamiltonianTransverse(oddSpins, bonds, N, J, h);
        eigtemp2=eigs(Symmetric(HtempOdd), nev=2^(N*N-1));
        HtempOdd=nothing;
        oddSpins=nothing;
        GC.gc();
        append!(eigenvalues, eigtemp[1]);
        append!(eigenvalues, eigtemp2[1]);
        append!(eigenvectors, eigtemp[2]);
        append!(eigenvectors, eigtemp2[2]);

    end

    append!(eigensystem, eigenvalues);
    append!(eigensystem, eigenvectors);

    return eigensystem;

end

#refStates, N, psector, bonds, refStatesMap
function calculateEigensystemHeisenbergMomentum(N, J, bonds)
    println("TIMING GENERATION REF STATES");
    @time begin
    refStatesData=referenceStates(N);
end
    println("TE SIZEH", length(refStatesData[1]));
    #list of all viable reference states
    refStates=refStatesData[1];
    #THIS maps a number to its state object, which contains info about ref state and shit
    refStatesMap=refStatesData[2];
    eigenvalues=Any[];
    for i=0:N-1
        println("CURR MOMENTUM ", i);
        viableSt=findViableRefStates(i, N, refStates);
        #viable st contains info about the numbering of each viable state
        Htemp::SparseMatrixCSC{Complex{Float64}}=constructHamiltonianHeisenbergMomentum(viableSt, N, i, bonds, refStatesMap);
        println(Htemp);
        println(Htemp[1,1]);
        eigtemp=eigs(Htemp, nev=100);
        append!(eigenvalues,eigtemp);
    end
    return eigenvalues;
end




#refStates, N, psector, bonds, refStatesMap
function calculateEigensystemHeisenbergMomentum2d(N, J, bonds)
    eigensystem=Any[];
    eigenvectors=Any[];
    count=0;
    println("TIMING GENERATION REF STATES");
    @time begin
    refStatesData=referenceStatesXY(N);
    end
    println("TE SIZEH", length(refStatesData[1]));
    #list of all viable reference states
    refStates=refStatesData[1];
    #THIS maps a number to its state object, which contains info about ref state and shit
    refStatesMap=refStatesData[2];
    momenta::Array{momentum}=generateAllMomenta(N);
    eigenvalues=Any[];
    for i=1:length(momenta)
        println("curr momentum: ", momenta[i]);
        local pt::momentum=momenta[i];
        viableSt=getViableStates2d(pt.px, pt.py, N, refStates);
        count+=length(viableSt[1]);
        #viable st contains info about the numbering of each viable state
        Htemp::SparseMatrixCSC{Complex{Float64}}=constructHamiltonianHeisenbergMomentum2d(viableSt, N, momenta[i], bonds, refStatesMap);
        println(Htemp);
        eigtemp=eigs(Htemp);
        append!(eigenvalues,eigtemp[1]);
        append!(eigenvectors, eigtemp[2])
    end
    push!(eigensystem, eigenvalues);
    push!(eigensystem, eigenvectors);
    println("TOTALTAOTA", count);
    return eigensystem;
end



#refStates, N, psector, bonds, refStatesMap
function calculateEigensystemHeisenbergReflection(N, J, bonds)
    eigensystem=Any[];
    eigenvectors=Any[];
    reflections=generateAllReflections();
    refStates=getReferenceStatesReflection(N);
    refStatesMap=refStates[2];
    viableStatesPos=findViableRefStatesReflection(refStates[1], 1, N);
    viableStatesNeg=findViableRefStatesReflection(refStates[1], -1, N);
#refStates, N, momentum, bonds, refStatesMap
    HPos=constructHamiltonianHeisenbergReflection(viableStatesPos,N, reflections,bonds, refStatesMap);
    HNeg=constructHamiltonianHeisenbergReflection(viableStatesNeg,N, reflections,bonds, refStatesMap);
    push!(eigensystem, eigenvalues);
    push!(eigensystem, eigenvectors);
    return eigensystem;
end
