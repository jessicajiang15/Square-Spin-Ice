include("hamiltonianGeneration.jl")
using LinearAlgebra
using KrylovKit
using Arpack

function calculateEigensystemHeisenberg(N, J, bonds)
    @time begin
        eigensystem::Array{Any}=Any[];

    eigenvalues::Array{Any}=Any[];
    eigenvectors::Array{Any}=Any[];
    for i=0:N*N
        spinUps::Array{Any}=singleOutNUpSpins(i, 2^(N*N));
        Htemp::SparseMatrixCSC{Float64}=constructHamiltonianHamiltonian(spinUps, bonds, N, J);
        #println(Htemp);
        if i==0||i==N*N
            append!(eigenvalues, Htemp[1, 1]);
            continue;
        end
        eigtemp=eigs(Htemp, nev=16);
        append!(eigenvalues, eigtemp[1]);
        append!(eigenvectors, eigtemp[2]);
    end
end
push!(eigensystem, eigenvalues);
push!(eigensystem, eigenvectors);

return eigensystem;
end



function calculateEigensystemTransverse(N, J, h, bonds)
    @time begin
        eigensystem::Array{Any}=Any[];
    eigenvalues::Array{Any}=Any[];
    eigenvectors::Array{Any}=Any[];
    evenSpins::Array{Any}=singleOutEvenOddSpins(true, 2^(N*N), N);
        oddSpins::Array{Any}=singleOutEvenOddSpins(false, 2^(N*N), N);
        println("STARTING EVEN");
        HtempEven::SparseMatrixCSC{Float64}=constructHamiltonianTransverse(evenSpins, bonds, N, J,h);

        println("STARTING ODD");
        HtempOdd::SparseMatrixCSC{Float64}=constructHamiltonianTransverse(oddSpins, bonds, N, J, h);
        println(HtempEven);
        println(HtempOdd);

        eigtemp=eigs(HtempEven);
        eigtemp2=eigs(HtempOdd);
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
    refStates=refStatesData[1];
    #THIS maps a number to its state object, which contains info about ref state and shit
    refStatesMap=refStatesData[2];
    eigenvalues=Any[];
    for i=-N*N/2+1:N/2
        println("CURR MOMENTUM ", i);
        viableSt=findViableRefStates(i, N, refStates);
        #viable st contains info about the numbering of each viable state
        Htemp::SparseMatrixCSC{Float64}=constructHamiltonianHeisenbergMomentum(viableSt, N, i, bonds, refStatesMap);
        eigtemp=eigs(Htemp);
        append!(eigenvalues,eigtemp);
    end
    return eigenvalues;
end
