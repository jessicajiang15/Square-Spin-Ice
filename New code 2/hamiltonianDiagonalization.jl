include("hamiltonianGeneration.jl")
using LinearAlgebra
using KrylovKit
using Arpack

function calculateEigensystemHeisenberg(N, J, bonds, eigmethod, num)
        eigensystem::Array{Any}=Any[];

    eigenvalues::Array{Any}=Any[];
    eigenvectors::Array{Any}=Any[];
    for i=0:N*N
        spinUps::Array{Any}=singleOutNUpSpins(i, 2^(N*N));
        @time begin
        #println("next:", i);
        Htemp=constructHeisenbergHamiltonian(spinUps, bonds, N*N, J, eigmethod);
    end
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
        n::Int=0;
        if(num=="all")
            n=length(spinUps[1]);
        elseif(num=="one")
            n=1;
        else
            n=16;
        end
        if(eigmethod=="full")
            println("yes", typeof(Htemp));
            eigtemp=eigen(Hermitian(Htemp));
            append!(eigenvalues, eigtemp.values);
            append!(eigenvectors, eigtemp.vectors);
        elseif(eigmethod=="lanczos")
            values, vecs, info=eigsolve(Htemp, n, :SR; krylovdim=100, ishermitian=true);
            append!(eigenvalues, values);
            append!(eigenvectors, vecs);
        else
            println(typeof(Htemp));
            println(n);
            eigtemp=eigs(Symmetric(Htemp), nev=n);
            #eigtemp=eigen(Htemp);
            #println(eigtemp);
            println(length(eigtemp[1]));
            append!(eigenvalues, eigtemp[1]);
            append!(eigenvectors, eigtemp[2]);

            #append!(eigenvalues, eigtemp.values);
            #append!(eigenvectors, eigtemp.vectors);
        end


        end

        push!(eigensystem, eigenvalues);
        push!(eigensystem, eigenvectors);
        return eigensystem;

end

#let N b e how many sites??
function calculateEigensystemTransverse(N, J, J2, h, bonds,eigmethod, num, hbar, width)
    #println("timing matrix generation");
    @time begin
    randomList=generateRandomh(hbar, width, bonds);
    #println(randomList);
    #println(randomList);
    #the info contains the states
    theInfo::Array{Any}=Any[];
    #contains a map of the states
    theInfo2::Array{Any}=Any[];
    eigensystem::Array{Any}=Any[];
    eigenvalues::Array{Any}=Any[];
    eigenvectors::Array{Any}=Any[];
    evenSpins=singleOutEvenOddSpins(true, 2^(N));
    #println("length", length(evenSpins[1]));

    push!(theInfo, evenSpins[1]);
    push!(theInfo2, evenSpins[2]);
    oddSpins=singleOutEvenOddSpins(false, 2^(N));
        #println("STARTING EVEN");
        @time begin
        HtempEven=constructTransverseHamiltonian(evenSpins, bonds, N, J, J2, eigmethod, randomList);
        #println("H[1,1]", HtempEven[1, 1]);
    end
    #println("size even", size(HtempEven));

    n::Int=0;
    if(num=="all")
        n=length(evenSpins[1]);
    elseif(num=="one")
        n=1;
    else
        n=16;
    end

        if(eigmethod=="full")
            eigtemp=eigen(Hermitian(HtempEven));
            append!(eigenvalues, eigtemp.values);
            append!(eigenvectors, eigtemp.vectors);
            #println("finished even eigenvalues");
        elseif(eigmethod=="lanczos")
            local vals=nothing
            local vecs=nothing
            local info=nothing;
            while(vals==nothing)
                try
                    vals, vecs, info=eigsolve(HtempEven, 1, :SR; krylovdim=200, ishermitian=true);
                catch
                    println("failed")
                end
            end
            println("even eigenvalues: ", vals);
            push!(eigenvalues, vals[1]);
            push!(eigenvectors, vecs[1]);
        else
            eigtemp=eigs(HtempEven, nev=n);
            #eigtemp=eigen(Htemp);
            #println(eigtemp);
            #println(length(eigtemp[1]));
            append!(eigenvalues, eigtemp[1]);
            append!(eigenvectors, eigtemp[2]);
            #append!(eigenvalues, eigtemp.values);
            #append!(eigenvectors, eigtemp.vectors);
        end
        #eigtemp=eigs(Symmetric(HtempEven), nev=2^(N*N-1));

        push!(theInfo, oddSpins[1]);
        push!(theInfo2, oddSpins[2]);

        #println("STARTING ODD");
        @time begin
        HtempOdd=constructTransverseHamiltonian(oddSpins, bonds, N, J, J2, eigmethod, randomList);
    end


    #println("size odd", size(HtempOdd));


        if(eigmethod=="full")
            eigtemp=eigen(Hermitian(HtempOdd));
            append!(eigenvalues, eigtemp.values);
            append!(eigenvectors, eigtemp.vectors);
        #    println("finished odd eigenvalues");
        elseif(eigmethod=="lanczos")
            local vals=nothing
            local vecs=nothing
            local info=nothing;
            while(vals==nothing)
                try
                    vals, vecs, info=eigsolve(HtempOdd, 1, :SR; krylovdim=200, ishermitian=true);
                catch
                    println("failed")
                end
            end
            println("odd eigenvalues: ", vals);

            push!(eigenvalues, vals[1]);
            push!(eigenvectors, vecs[1]);
            #println("vecs", vecs[1]);
            #println("length vecs, ", length(vecs[1]));
        else
            eigtemp=eigs((HtempOdd), nev=n);
            #eigtemp=eigen(Htemp);
            #println(eigtemp);
            #println(length(eigtemp[1]));
            append!(eigenvalues, eigtemp[1]);
            append!(eigenvectors, eigtemp[2]);
            #append!(eigenvalues, eigtemp.values);
            #append!(eigenvectors, eigtemp.vectors);
        end
    end

    push!(eigensystem, eigenvalues);
    push!(eigensystem, eigenvectors);
    push!(eigensystem, theInfo);
    push!(eigensystem, theInfo2);
    return eigensystem;
end



#refStates, N, psector, bonds, refStatesMap
function calculateEigensystemHeisenbergMomentum(N, J, bonds, eigmethod, num)
    println("TIMING GENERATION REF STATES");
    @time begin
    refStatesData=referenceStates(N);
end
    eigensystem=Any[];
    println("TE SIZEH", length(refStatesData[1]));
    #list of all viable reference states
    refStates=refStatesData[1];
    #THIS maps a number to its state object, which contains info about ref state and shit
    refStatesMap=refStatesData[2];
    eigenvalues=Any[];
    eigenvectors=Any[];
    for i=0:N-1
        println("CURR MOMENTUM ", i);
        viableSt=findViableRefStates(i, N, refStates);
        #viable st contains info about the numbering of each viable state
        Htemp=constructHamiltonianHeisenbergMomentum(viableSt, N*N, i, bonds, refStatesMap, eigmethod);
        println(Htemp);
        println(Htemp[1,1]);
        n::Int=0;
        if(num=="all")
            n=length(spinUps[1]);
        elseif(num=="one")
            n=1;
        else
            n=16;
        end
        if(eigmethod=="full")
            eigtemp=eigen(Hermitian(Htemp));
            append!(eigenvalues, eigtemp.values);
            append!(eigenvectors, eigtemp.vectors);
        elseif(eigmethod=="lanczos")
            values, vecs, info=eigsolve(Htemp, n, :SR; krylovdim=100, ishermitian=true);
            append!(eigenvalues, values);
            append!(eigenvectors, vecs);
        else
            eigtemp=eigs((Htemp), nev=n);
            #eigtemp=eigen(Htemp);
            #println(eigtemp);
            println(length(eigtemp[1]));
            println(eigtemp[1]);
            append!(eigenvalues, eigtemp[1]);
            append!(eigenvectors, eigtemp[2]);

            #append!(eigenvalues, eigtemp.values);
            #append!(eigenvectors, eigtemp.vectors);
        end
    end
    push!(eigensystem, eigenvalues);
    push!(eigensystem, eigenvectors);
    return eigensystem;
end




#refStates, N, psector, bonds, refStatesMap
function calculateEigensystemHeisenbergMomentum2d(N, J, bonds, eigmethod, num)
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
        Htemp=constructHamiltonianHeisenbergMomentum2d(viableSt, N*N, momenta[i], bonds, refStatesMap, eigmethod);
        println(Htemp);
        println("IS HERMITIAN: ", isHermitian(Htemp));
        n::Int=0;
        if(num=="all")
            n=length(viableSt[1]);
        elseif(num=="one")
            n=1;
        else
            n=16;
        end
        if(eigmethod=="full")
            eigtemp=eigen(Hermitian(Htemp));

            append!(eigenvalues, eigtemp.values);
            println("THE MOMENTA", momenta[i].px, ", ", momenta[i].py)
            println(eigtemp.values);
            append!(eigenvectors, eigtemp.vectors);
        elseif(eigmethod=="lanczos")
            values, vecs, info=eigsolve(Htemp, n, :SR; krylovdim=100, ishermitian=true);
            append!(eigenvalues, values);
            append!(eigenvectors, vecs);
        else
            eigtemp=eigs((Htemp), nev=n);
            #eigtemp=eigen(Htemp);
            #println(eigtemp);
            println(length(eigtemp[1]));
            #println("eigenvalues of ", momenta[i].px, ", ", momenta[1].py);
            #println(eigtemp[1]);
            append!(eigenvalues, eigtemp[1]);
            append!(eigenvectors, eigtemp[2]);

            #append!(eigenvalues, eigtemp.values);
            #append!(eigenvectors, eigtemp.vectors);
        end
    end
    length(eigenvalues);
    push!(eigensystem, eigenvalues);
    push!(eigensystem, eigenvectors);
    println("TOTALTAOTA", count);
    return eigensystem;
end



#refStates, N, psector, bonds, refStatesMap
function calculateEigensystemHeisenbergReflection(N, J, bonds, eigmethod, num)
    eigensystem=Any[];
    eigenvectors=Any[];
    eigenvalues=Any[];
    reflections=generateAllReflections();
    refs=getReferenceStatesReflection(N);
    lambdas=generateAllReflections();
    println("how many refs", length(refs[1]));
    total=0;
    refStatesMap=refs[2];
    for i=1:length(reflections)
        viableStates=findViableRefStatesReflection(refs[1], reflections[i], N);
        println("THE LENGTH",length(viableStates[1]));
        total+=length(viableStates[1]);
        n::Int=0;
        if(num=="all")
            n=length(viableStates[1]);
        elseif(num=="one")
            n=1;
        else
            n=100;
        end
        Htemp=constructHamiltonianHeisenbergReflection(viableStates,N*N, reflections[i],bonds, refStatesMap, eigmethod);
        println("is hermitian: ", isHermitian(Htemp));
        if(eigmethod=="full")
            eigtemp=eigen(Htemp);
            append!(eigenvalues, eigtemp.values);
            append!(eigenvectors, eigtemp.vectors);
        elseif(eigmethod=="lanczos")
            values, vecs, info=eigsolve(Htemp, n, :SR; krylovdim=100, ishermitian=true);
            append!(eigenvalues, values);
            append!(eigenvectors, vecs);
        else
            eigtemp=eigs((Htemp), nev=n);
            #eigtemp=eigen(Htemp);
            #println(eigtemp);
            println(length(eigtemp[1]));
            append!(eigenvalues, eigtemp[1]);
            append!(eigenvectors, eigtemp[2]);
            #append!(eigenvalues, eigtemp.values);
            #append!(eigenvectors, eigtemp.vectors);
        end
    end
    println("TOTAOKL", total);
    push!(eigensystem, eigenvalues);
    push!(eigensystem, eigenvectors);
    return eigensystem;
end





#let N b e how many sites??
function calculateEigensystemSusceptibility(N, J, J2, h, h2, bonds,eigmethod, num, hbar, hbar2, width, sites2)
    #println("timing matrix generation");
    @time begin
    randomList=generateRandomh(hbar, width, bonds);
    randomList2=generateRandomh(hbar2, width, bonds);
    theInfo::Array{Any}=Any[];
    eigensystem::Array{Any}=Any[];
    eigenvalues::Array{Any}=Any[];
    eigenvectors::Array{Any}=Any[];
    local Htemp;
    push!(theInfo, 0:2^(N)-1);
    #println("timing h");
    Htemp=constructSusceptibilityHamiltonian(bonds, N, J, J2, eigmethod, randomList, randomList2, sites2);

    n::Int=0;
    if(num=="all")
        n=length(evenSpins[1]);
    elseif(num=="one")
        n=1;
    else
        n=16;
    end
    if(eigmethod=="full")
        eigtemp=eigen(Hermitian(Htemp));
        append!(eigenvalues, eigtemp.values);
        append!(eigenvectors, eigtemp.vectors);
        #println("finished odd eigenvalues");
    elseif(eigmethod=="lanczos")
        local values=nothing
        local vecs=nothing
        local info=nothing;
        while(values==nothing)
            try
                values, vecs, info=eigsolve(Htemp, 1, :SR; krylovdim=200, ishermitian=true, tol=10^(-16));
            catch
                println("failed")
            end
        end
        #println(values);
        #println(values);
        append!(eigenvalues, values[1]);
        append!(eigenvectors, vecs[1]);
    else
        eigtemp=eigs((Htemp), nev=n);
        #eigtemp=eigen(Htemp);
        #println(eigtemp);
        #println(length(eigtemp[1]));
        append!(eigenvalues, eigtemp[1]);
        append!(eigenvectors, eigtemp[2]);
        #append!(eigenvalues, eigtemp.values);
        #append!(eigenvectors, eigtemp.vectors);
    end

push!(eigensystem, eigenvalues);
push!(eigensystem, eigenvectors);
push!(eigensystem, theInfo);
end
return eigensystem;
end



#numSites, J1, J2, hs, bonds,"lanczos", "one"
#takes in a list of fields at each site and applies that field
function calculateEigensystemTransverse(N, J, J2, hs::Vector{Float64}, bonds, eigmethod, num)
    #println("timing matrix generation");
    @time begin
    #println(randomList);
    #println(randomList);
    #the info contains the states
    theInfo::Array{Any}=Any[];
    #contains a map of the states
    theInfo2::Array{Any}=Any[];
    eigensystem::Array{Any}=Any[];
    eigenvalues::Array{Any}=Any[];
    eigenvectors::Array{Any}=Any[];
    evenSpins=singleOutEvenOddSpins(true, 2^(N));
    #println("length", length(evenSpins[1]));

    push!(theInfo, evenSpins[1]);
    push!(theInfo2, evenSpins[2]);
    oddSpins=singleOutEvenOddSpins(false, 2^(N));
        #println("STARTING EVEN");
        @time begin
        HtempEven=constructTransverseHamiltonian(evenSpins, bonds, N, J, J2, eigmethod, hs);
        #println("H[1,1]", HtempEven[1, 1]);
    end
    #println("size even", size(HtempEven));

    n::Int=0;
    if(num=="all")
        n=length(evenSpins[1]);
    elseif(num=="one")
        n=1;
    else
        n=16;
    end

        if(eigmethod=="full")
            eigtemp=eigen(Hermitian(HtempEven));
            append!(eigenvalues, eigtemp.values);
            append!(eigenvectors, eigtemp.vectors);
            #println("finished even eigenvalues");
        elseif(eigmethod=="lanczos")
            local vals=nothing
            local vecs=nothing
            local info=nothing;
            while(vals==nothing)
                try
                    vals, vecs, info=eigsolve(HtempEven, 1, :SR; krylovdim=200, ishermitian=true);
                catch
                    println("failed")
                end
            end
            println("even eigenvalues: ", vals);
            push!(eigenvalues, vals[1]);
            push!(eigenvectors, vecs[1]);
        else
            eigtemp=eigs(HtempEven, nev=n);
            #eigtemp=eigen(Htemp);
            #println(eigtemp);
            #println(length(eigtemp[1]));
            append!(eigenvalues, eigtemp[1]);
            append!(eigenvectors, eigtemp[2]);
            #append!(eigenvalues, eigtemp.values);
            #append!(eigenvectors, eigtemp.vectors);
        end
        #eigtemp=eigs(Symmetric(HtempEven), nev=2^(N*N-1));

        push!(theInfo, oddSpins[1]);
        push!(theInfo2, oddSpins[2]);

        #println("STARTING ODD");
        @time begin
        HtempOdd=constructTransverseHamiltonian(oddSpins, bonds, N, J, J2, eigmethod, hs);
    end


    #println("size odd", size(HtempOdd));


        if(eigmethod=="full")
            eigtemp=eigen(Hermitian(HtempOdd));
            append!(eigenvalues, eigtemp.values);
            append!(eigenvectors, eigtemp.vectors);
        #    println("finished odd eigenvalues");
        elseif(eigmethod=="lanczos")
            local vals=nothing
            local vecs=nothing
            local info=nothing;
            while(vals==nothing)
                try
                    vals, vecs, info=eigsolve(HtempOdd, 1, :SR; krylovdim=200, ishermitian=true);
                catch
                    println("failed")
                end
            end
            println("odd eigenvalues: ", vals);

            push!(eigenvalues, vals[1]);
            push!(eigenvectors, vecs[1]);
            #println("vecs", vecs[1]);
            #println("length vecs, ", length(vecs[1]));
        else
            eigtemp=eigs((HtempOdd), nev=n);
            #eigtemp=eigen(Htemp);
            #println(eigtemp);
            #println(length(eigtemp[1]));
            append!(eigenvalues, eigtemp[1]);
            append!(eigenvectors, eigtemp[2]);
            #append!(eigenvalues, eigtemp.values);
            #append!(eigenvectors, eigtemp.vectors);
        end
    end

    push!(eigensystem, eigenvalues);
    push!(eigensystem, eigenvectors);
    push!(eigensystem, theInfo);
    push!(eigensystem, theInfo2);
    return eigensystem;
end





#to get first eigenvalue: use temp[1][1]
#to get first eigenvector: use temp[2][1]
#to get states list: use temp[3][1]
function calculateEigensystemTransverseNoSymmetry(N, J, J2, hs::Vector{Float64},hs2::Vector{Float64}, bonds,eigmethod, num, h1orh2)
    #println(randomList);
    #the info contains the states
    theInfo::Array{Any}=Any[];
    eigensystem::Array{Any}=Any[];
    eigenvalues::Array{Any}=Any[];
    eigenvectors::Array{Any}=Any[];
    local Htemp;
    push!(theInfo, 0:2^(N)-1);
    println("timing h");
    @time begin
        if(h1orh2=="H1")
            #H1= field is in x direction
            Htemp=constructTransverseHamiltonianNoSymmetrySxMeanField(bonds, N, J, J2, eigmethod, hs, hs2);
        else
            #H2= field is in the z direction
            Htemp=constructTransverseHamiltonianNoSymmetrySz(bonds, N, J, J2, eigmethod, hs);
        end
    end
    n::Int=0;
    if(num=="all")
        n=2^N;
    elseif(num=="one")
        n=1;
    else
        n=5;
    end
    if(eigmethod=="full")
        eigtemp=eigen(Hermitian(Htemp));
        append!(eigenvalues, eigtemp.values);
        append!(eigenvectors, eigtemp.vectors);
        #println("finished odd eigenvalues");
    elseif(eigmethod=="lanczos")
        println("starting lanczos")
        local values=nothing
        local vecs=nothing
        local info=nothing;
        while(values==nothing)
            try
                values, vecs, info=eigsolve(Htemp, 1, :SR; krylovdim=200, ishermitian=true, tol=10^(-16));
            catch
                println("lanczos failed, retry")
            end
        end
        append!(eigenvalues, values[1]);
        append!(eigenvectors, vecs[1]);
    else
        eigtemp=eigs((Htemp), nev=n, which="SR");
        append!(eigenvalues, eigtemp[1]);
        append!(eigenvectors, eigtemp[2]);
    end

push!(eigensystem, eigenvalues);
push!(eigensystem, eigenvectors);
push!(eigensystem, theInfo);
return eigensystem;

end

