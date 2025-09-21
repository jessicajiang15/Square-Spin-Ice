using Plots
include("necessaryStructures.jl")
using JLD
using .Threads

function readFromFile(filename)
    list::Array{Complex{Float64}}=Complex{Float64}[];
    f=open(filename*".txt", "r");
    while !eof(f)
        s=readLine(f);
        push!(list, s);
    end
    return list;
end

function partition(beta, E)
    sum=0;
    for i=1:length(E)
        sum+=exp(E[i]*beta);
    end
    return sum;
end

function calculateExpE(beta, E)
    sum=0;
    for i=1:length(E)
        sum+=E[i]*exp(E[i]*beta);
    end
    sum/=partition(beta, E);
    return sum;
end

function calculateExpE(beta, E, partition)
    sum=0;
    for i=1:length(E)
        sum+=E[i]*exp(E[i]*beta);
    end
    sum/=partition;
    return sum;
end

function calcualteExpE(beta, E, partition)
    sum=0;
    for i=1:length(E)
        sum+=(E[i]^2)*exp(E[i]*beta);
    end
    sum/=partition;
    return sum;
end


function calculateExpESq(beta, E)
    sum=0;
    for i=1:length(E)
        sum+=(E[i]^2)*exp(E[i]*beta);
    end
    sum/=partition(beta, E);
    return sum;
end

function calculateExpESq(beta, E, partition)
    sum=0;
    for i=1:length(E)
        sum+=(E[i]^2)*exp(E[i]*beta);
    end
    sum/=partition;
    return sum;
end



function calculateHeatCapacity(beta, E, E2)
    return (beta^2)*(E2-E^2);
end

function calculateG(beta, partition)
    return (-1/beta)*log(partition);
end

function calculateEntropy(beta, G, E)
    return beta*(G-E);
end

function sumDigits(num, N)
    n=num;
    n = (n & 0x55555555) + (n >>  1 & 0x55555555);
    n = (n & 0x33333333) + (n >>  2 & 0x33333333);
    n = (n & 0x0F0F0F0F) + (n >>  4 & 0x0F0F0F0F);
    n = (n & 0x00FF00FF) + (n >>  8 & 0x00FF00FF);
    n = (n & 0x0000FFFF) + (n >> 16 & 0x0000FFFF);
return n-N*N+n;
end

function getAllRelevantQuantities(bmin, bmax, bstep, eigenvalues)
    all=Any[];
    heatCapacities=Float64[];
    gs=Float64[];
    expESqs=Float64[];
    expEs=Float64[];
    partitions=Float64[];
    entropies=Float64[];
    betas=Float64[];
    b=bmin;

    while b<=bmax
        part=partition(b, eigenvalues);
        push!(partitions, part);
        expE=calculateExpE(b, eigenvalues, part);
        push!(expEs, expE);
        expESq=calculateExpESq(b, eigenvalues, part);
        push!(expESqs, expESq);
        heatCap=calculateHeatCapacity(b, expE, expESq);
        push!(heatCapacities, heatCap);
        g=calculateG(b, part);
        push!(heatCapacities, heatCap);
        entropy=calculateEntropy(b, g, expE);
        push!(entropies, entropy);
        b+=bstep;
    end
    push!(all, betas);
    push!(all, partitions);
    push!(all, expEs);
    push!(all, expESqs);
    push!(all, heatCapacities);
    push!(all, gs);
    push!(all, entropies);
    return all;
end

function calculateMagnetization(N)
    sum=0;
    for i=1:2^(N*N)
        #fill this in later
    end

end

#=
push!(all, betas);
push!(all, partitions);
push!(all, expEs);
push!(all, expESqs);
push!(all, heatCapcities);
push!(all, gs);
push!(all, entropies);
heat cap=5
g=6
entropy=7
=#
function theplot(list, i::Int, str)
    println("IM PLOTTING")
    plot(list[1],list[i])
    #
    png("Users/Jessica/git/Square Spin Ice/plot.png")
end

function plotH(list, list2, str)
    println("IM PLOTTING")
    plot(list,list2, title = string(str, " vs. h"))
    savefig("Users/Jessica/git/Square Spin Ice/plot.png")
end

function plotQuantity(bmin, bmax, bstep, eigenvalues, i)
    all=getAllRelevantQuantities(bmin, bmax, bstep, eigenvalues);
    str::String="";
    if(i==5)
        str="Heat Capacity";
    elseif(i==6)
        str="Free Energy";
    else
        str="Entropy";
    end
    theplot(all, i, str);
end


function plotQuantity(all, i)
    str::String="";
    if(i==5)
        str="Heat Capacity";
    elseif(i==6)
        str="Free Energy";
    else
        str="Entropy";
    end
    theplot(all, i, str);
end

#N IS TOTAL NUMBER OF STATES
function calculateSz(eigenvector, N)
    sum=0;
    for i=0:2^(N)-1
        temp=countBits(i);
        sum+=(temp-(N-temp))*abs2(eigenvector[i+1]);
    end
    return (1/(N))*(1/2)*sum;
end
#N IS TOTAL NUMBER OF STATES
function calculateSz(eigenvector, states, N)
    sum=0;
    for i=1:length(states)
        temp=countBits(states[i]);
        sum+=(1/2)*(temp-(N-temp))*conj(eigenvector[i])*eigenvector[i];
    end
    return -(1/(N))*sum;
end
#unused
function szMatrix(eigenvector, states, N)
    m::Matrix{Float64}=zeros(length(eigenvector), length(eigenvector));
    for i=length(states)
        temp=countBits(states[i]);
        m[i, i]+=(1/2)*(temp-(N*N-temp));
    end
    return conj.(eigenvector')*m*eigenvector*(1/(N*N));
end
#N IS TOTAL NUMBER OF STATES
function calculateSx(eigenvector, N)
    sum=0;
    for i=0:2^(N)-1
        for j=1:N
            b=flipBit(j-1, i);
            sum+=eigenvector[i+1]*conj(eigenvector[b+1]);
        end
    end
    return (1/2)*sum;
end


function calculateSx(eigenvector, eigenvector2, N)
    sum=0;
    for i=0:2^(N)-1
        for i=1:N
            b=flipBit(i-1, i);
            sum+=eigenvector[i+1]*conj(eigenvector2[b+1]);
        end
    end
    return (1/2)*sum;
end

#let state be map from state -> position in array
#N IS TOTAL NUMBER OF STATES
function calculateSx(eigenvector, states, map, N)
    sum=0;
    for i=1:length(states)
        for j=1:N
            b=flipBit(j-1, states[i]);
            if(b in keys(map))
                c=map[b];
                sum+=eigenvector[i]*conj(eigenvector[c]);
            end
        end
    end
    return (1/2)*sum;
end


function isNeel(state, sector, N)
    #println("sector", sector);
    one=getTi(sector[1]-1, state);
    two=getTi(sector[2]-1, state);
    three=getTi(sector[3]-1, state);
    four=getTi(sector[4]-1, state);
    return one!=two&&three!=one&&four!=three&&four!=two;
end


function indiciesAllSquares(N)
    ind=Int[];
    for i=1:N*(N-1)
        if(i%4!=0)
            push!(ind, i);
        end
    end
    return ind;
end


#let state be map from state -> position in array
function calculateFlippabilityExp(eigenvector, states, N)
    sum=0;
    for i=1:length(states)
        temp=0;
        for j=1:N*(N-1)
            if(j%4!=0&&isNeel(states[i], j, N))
                temp+=1;
            end
        end
        sum+=temp*abs2(eigenvector[i]);
    end
    return sum;
end

function calculateFlippabilityExp(eigenvector, states, sector, N)
    sum=0;
    for i=1:length(states)
        temp=0;
            if(isNeel(states[i], sector, N))
                sum+=temp*abs2(eigenvector[i]);
            end
    end
    return sum;
end

function neg(i::Int, j::Int)
    return (-1)^(i%2==0 ? 1 : -1 + j%2==0 ? 1 : -1);
end

function calculateStaggeredFlippability(eigenvector, states, squareIndicies, N)
    sum=0;
    #square indicies= list of list of indicies of the plaquettes, the first element in each list is the leading index
    for i=1:length(states)
        for j=1:length(squareIndicies)
            for z=j+1:length(squareIndicies)

                product=calculateFlippabilityExp(eigenvector, states, squareIndicies[j][1])*calculateFlippabilityExp(eigenvector, states, squareIndicies[z][1]);
                sum+=product*(-1)^(getPlaquetteNumber(squareIndicies[j][1], N)*getPlaquetteNumber(squareIndicies[z][1], N));
            end
        end
    end
    return sum*(N/2)^(-4);
end

function calculateSxMomentum(ref::refState2d, momentum::momentum, phases::Dict{Int, Complex{Float64}}, N)
    eigenvector::Vector{Complex{Float64}}=[];
    list::Array{Int}=Int[];
    Na=1/sqrt(ref.numUniqueSt*abs2(phases[ref.state]));
    push!(list, ref.state);
    for x=0:N-1
        for y=0:N-1
            if(x==0&&y==0)
                continue;
            end
            temp=rotateBits(x, y, state, N);
                push!(list, temp);
                push!(eigenvector, (Na)*exp(-im*2pi(momentum.px*x+momentum.py*y)/N));
        end
    end
    return calculateSx(eigenstate, list, N);

end



function calculateSxMomentumFull(eigenstate, refStates, momentum::momentum, phases::Dict{Int, Complex{Float64}}, N)
    sum=0;
    for i=1:length(eigenstate)
        for j=1:length(refStates)
            calculateSzMomentum(refStates[j], momentum, phases, N);
            sum+=abs2(eigenstate[i])*calculateSzMomentum(refStates[j], momentum, phases, N);
        end
    end
    return calculateSx(eigenstate, list, N);
end

function getLowestLyingStates(eigenvalues, eigenvectors)
    eigensystem=Any[];
    min=eigenvalues[1];
    index=1;
    #println("length eigenvalues, ", length(eigenvalues));
    for i=1:length(eigenvalues)
        if(eigenvalues[i]<min)
            min=eigenvalues[i];
            index=i;
        end
    end
    push!(eigensystem, min);
    push!(eigensystem, eigenvectors[index]);
    push!(eigensystem, index);
    return eigensystem;
end

function isIn(n, list)
     for i=1:length(list)
         if(n==list[i])
             return i;
         end
     end
     return -1;
end

#inner product <eigenvector|eigenvector2>, assumes states are separate
function innerProduct(eigenvector, states, eigenvector2, states2)
    sum=0;
    for i=1:length(states)
        index=findfirst(==(states[i]), states2)
        #binarySearchIt(states[i], states2);
        if(index!==nothing)
            sum+=conj(eigenvector[i])*eigenvector2[index];
        end
    end
    return sum;
end

#assumes the two states are the same
function innerProduct(eigenvector, eigenvector2, states)
    sum=0;
    for i=1:length(states)
           sum+=conj(eigenvector[i])*eigenvector2[i];
    end
    return sum;
end

function calculateFidelity(eigenvector, states, h, N, deltah, bonds, J, J2)
    temp=calculateEigensystemTransverse(N, J, J2, h+deltah, bonds,"lanczos", "one", h+deltah, 0);
    eigensystem=getLowestLyingStates(temp[1], temp[2]);
    eigenvector2=eigensystem[2];
    states2=temp[3][eigensystem[3]];
    #=
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    temp[3][eigensystem[3]]
    =#
    inner=abs(innerProduct(eigenvector, states, eigenvector2, states2));
    println("inner product: ", abs(innerProduct(eigenvector, states, eigenvector2, states2)));
    return 2*(1-inner)/(deltah^2);
end

#list is a list of all states

function generateFidelityList(hmin, hmax, num, N, J, J2, bonds)
    hstep=hmax/num-hmin/num;
    i=hmin;
    eigenvectors=Any[];
    states=Any[];
    all=Any[];
    println(J);

    while(i<=hmax+hstep)
        println("on i: ", i);
        temp=calculateEigensystemTransverse(N, J, J2, i, bonds,"lanczos", "one", i, 0);
        eigensystem=getLowestLyingStates(temp[1], temp[2]);
        eigenvector=eigensystem[2];
        state=temp[3][eigensystem[3]];
        push!(eigenvectors, eigenvector);
        push!(states, state);
        if(hstep<=0)
            break;
        end
        i+=hstep;
    end
    push!(all, eigenvectors);
    push!(all, states);
    return all;
end


function generateFidelityListH1(hmin, hmax, num, N, J, J2, bonds)
    hstep=hmax/num-hmin/num;
    i=hmin;
    eigenvectors=Any[];
    states=Any[];
    all=Any[];
    println(J);

    while(i<=hmax+hstep)
        println("on i: ", i);
        temp=calculateEigensystemTransverseNoSymmetry(N, J, J2, i, bonds,"lanczos", "one", i, 0, "H1");
        #eigensystem=getLowestLyingStates(temp[1], temp[2]);
        #eigenvector=eigensystem[2];
        push!(eigenvectors, temp[2]);
        push!(states, 0:2^(N)-1);
        if(hstep<=0)
            break;
        end
        i+=hstep;
    end
    push!(all, eigenvectors);
    push!(all, states);
    return all;
end

function calculateFidelity(hmin, hmax, num, N, J, J2, bonds)
    println("begin");
    println(J);
    fidelities=Any[];
    hstep=hmax/num-hmin/num;
    all=Any[];
    println("hlist begin");
    hlist=generateHListUniform(hmin, hmax, num);
    println("fidelity list begin");
    list=generateFidelityList(hmin, hmax, num, N, J, J2, bonds);
    println("fidelity list length", length(list));
    eigenvectors=list[1];
    states=list[2];
    for i=1:length(eigenvectors)-1
        inner=abs(innerProduct(eigenvectors[i], states[i], eigenvectors[i+1], states[i+1]));
        push!(fidelities,2*(1-inner)/(hstep^2));
    end
    return fidelities;
end




function calculateFidelity(hstep, stateList)
    fidelities=Any[];
    println("length", length(stateList));
    for i=1:length(stateList)-1
        println("kalsa");
        eigenvector=stateList[i][1];
        state=stateList[i][2];
        eigenvector2=stateList[i+1][1];
        state2=stateList[i+1][2];
        inner=abs(innerProduct(eigenvector, state, eigenvector2, state2));
        push!(fidelities,2*(1-inner)/(hstep^2));
    end
    return fidelities;
end

function calculateFidelity(hstep, fidelityList, i)
    eigenvectors=fidelityList[1];
    states=fidelityList[2];
    inner=abs(innerProduct(eigenvectors[i], states[i], eigenvectors[i+1], states[i+1]));
    return 2*(1-inner)/(hstep^2)
end


function getSzPlaquette(state, plaquetteIndicies)
    sum=0;
    println("plq",plaquetteIndicies);
    for i=1:length(plaquetteIndicies)
        sum+=getTi(plaquetteIndicies[i], state) == 1 ? 1 : -1;
    end
    return (1/2)*sum;
end

function calculateSxIndicies(eigenvector, states, map, N, indicies)
    sum=0;
    for i=1:length(states)
        for j=1:length(indicies)
            b=flipBit(indicies[j]-1, states[i]);
            if(b in keys(map))
                println("yayayayayayaya");
                c=map[b];
                sum+=eigenvector[i]*conj(eigenvector[c]);
            end
        end
    end
    return (1/2)*sum;
end



#calculation: sum (ij) sxisxj -1^(i+j)
function calculateSPi(eigenvector, states, N, map)
    plaquettes=generateListsofPlaquetteIndicies(N);
    totalsum=0;
    for i=1:length(states)
        println("state: ", i);
        for j=1:length(plaquettes)
            for z=j+1:length(plaquettes)
                #println("im done 1");
                jSx=calculateSxIndicies(eigenvector, states, map, N, plaquettes[j]);
                zSx=calculateSxIndicies(eigenvector, states, map, N, plaquettes[z]);
                #println("im done 2");
                no1=getPlaquetteNumber(plaquettes[j][1], N)
                no2=getPlaquetteNumber(plaquettes[z][1], N)
                #println("im done 3");
                #plauettes[j][1] gives the leading index of the plaquette
                totalsum+=jSx*zSx*abs2(eigenvector[i])*(-1)^(no1+no2);
                #println("im done 4");
            end
        end
    end
    return totalsum;
end
#sum(sx sx)


function calculateSPiSz(eigenvector, states, N)
    plaquettes=generateListsofPlaquetteIndicies(N);
    totalsum=0;
    for i=1:length(states)
        println("state: ", i);
        for j=1:length(plaquettes)
            for z=j+1:length(plaquettes)
                #println("im done 1");
                jSz=getSzPlaquette(states[i], plaquettes[j]);
                zSz=getSzPlaquette(states[i], plaquettes[z]);
                #println("im done 2");
                no1=getPlaquetteNumber(plaquettes[j][1], N)
                no2=getPlaquetteNumber(plaquettes[z][1], N)
                #println("im done 3");
                #plauettes[j][1] gives the leading index of the plaquette
                totalsum+=jSz*zSz*abs2(eigenvector[i])*(-1)^(no1+no2);
                #println("im done 4");
            end
        end
    end
    return totalsum;
end



function calculateSPiSzNew(eigenvector, states, N)
    totalsum=0;
        for j=1:N
            for z=1:N
                for i=1:length(states)
                    jSz=getTi(j-1, states[i])==1 ? 1/2 : -1/2;
                    zSz=getTi(z-1, states[i])==1 ? 1/2 : -1/2;
                    totalsum+=jSz*zSz*abs2(eigenvector[i])*(-1)^(getPlaquetteNumber(z, round(sqrt(N)))+getPlaquetteNumber(j, round(sqrt(N))));
            end
        end
    end
    return totalsum/(N);
end





function calculateSPiSzNew(eigenvector, states, N, indicies)
    totalsum=0;
    println("length,", length(eigenvector))
    println("length, ", length(indicies));

    for i=1:length(states)
        #println("state: ", i);
        for j=1:N
            for z=1:N
                #println("im done 1");
                jSz=getTi(j-1, states[i])==1 ? 1/2 : -1/2;
                zSz=getTi(z-1, states[i])==1 ? 1/2 : -1/2;
                #println("im done 2");
                #println("im done 3");

                #plauettes[j][1] gives the leading index of the plaquette
                totalsum+=jSz*zSz*abs2(eigenvector[i])*(-1)^(indicies[z]+indicies[j]);
                #*(-1)^(i+j)
                #println("im done 4");
            end
        end
    end
    return totalsum/(N);
end



function calculateSPiSzNewAbs(eigenvector, states, N)
    totalsum=0;
        #println("state: ", i);
        for j=1:N
            for z=1:N
                tempsum=0;
                for i=1:length(states)
                #println("im done 1");
                    jSz=getTi(j-1, states[i])==1 ? 1/2 : -1/2;
                    zSz=getTi(z-1, states[i])==1 ? 1/2 : -1/2;
                    #println("im done 2");
                    #println("im done 3");
                    #plauettes[j][1] gives the leading index of the plaquette
                    tempsum+=(jSz*zSz)*(abs2(eigenvector[i]));
                    #*(-1)^(i+j)
                    #println("im done 4");
                end
                totalsum+=abs(tempsum);
            end
        end
    return totalsum/(N);
end


function calculateFlippabilityNew(eigenvector, states, N, plaquettes)
    #println("plaquettes", plaquettes);
    total=0;
    for i=1:length(states)
        for j=1:length(plaquettes)
            for z=1:length(plaquettes)
                #println(plaquettes[j]);
                if(isNeel(states[i], plaquettes[j], N)&&isNeel(states[i], plaquettes[z], N))
                    total+=abs2(eigenvector[i])*(-1^(j+z));
                end
            end
        end
    end
    return -total/(length(plaquettes));
end

#=
function generateHListUniform(J, J2, num)
    list::Array{Float64}=Float64[];
    i=(1/10)*J;
    interval=(J-i)/num;

    while(i<=J)
        push!(list, i);
        i+=interval;
    end
    return list;
end
=#

function generateHListUniformHalf(J, J2, num)
    list::Array{Float64}=Float64[];
    i=(1/10)*J;
    interval=(0.5*J-i)/num;

    while(i<=0.5*J)
        push!(list, i);
        i+=interval;
    end
    return list;
end

function generateHListLog(J, J2, num)
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


function generateHListUniform(hmin, hmax, num)
    inc=hmax/num-hmin/num;
    i=hmin;
    hlist=Float64[];
    while(i<=hmax)
        push!(hlist, i);
        i+=inc;
        if(inc==0)
            break;
        end
        println("inc", inc);
    end
    return hlist;
end



function generateHListUniformIncludeOne(hmin, hmax, num)
    inc=hmax/num-hmin/num;
    i=hmin;
    yes=false;
    hlist=Float64[];
    while(i<=hmax)
        push!(hlist, i);
        i+=inc;
        if(i>1&&i!=1&&!yes)
            push!(hlist, 1);
            yes=true;
        end
        if(inc==0)
            break;
        end
    end
    return hlist;
end



function storeInformation(bonds, N)
    plaquettes=generateListsofPlaquetteIndiciesFlip(N);

end


function calculateSusceptibility(E0, Eh, hs)
    return (E0-Eh)/(hs^2)
end

function calculateNeelOrderSz(eigenvector, states, N)
    totalsum=0;
        for j=1:N
                for i=1:length(states)
                    jSz=getTi(j-1, states[i])==1 ? 1/2 : -1/2;
                    totalsum+=jSz*abs2(eigenvector[i])*(-1)^(getPlaquetteNumber(j, round(sqrt(N))));
        end
    end
    return totalsum/(N);
end




function calculateNeelOrderSz(eigenvector, states, N, indicies)
    totalsum=0;
    for i=1:length(states)
        #println("state: ", i);
        for j=1:N
                #println("im done 1");
                jSz=getTi(j-1, states[i])==1 ? 1/2 : -1/2;
                totalsum+=jSz*abs2(eigenvector[i])*(-1)^(indicies[j]);
            end
    end
    return totalsum/(N);
end

function calculateFidelityPlain(eigenvector, states, h, N, deltah, bonds, J, J2)
    temp=calculateEigensystemTransverse(N, J, J2, h+deltah, bonds,"lanczos", "one", h+deltah, 0);
    eigensystem=getLowestLyingStates(temp[1], temp[2]);
    eigenvector2=eigensystem[2];
    states2=temp[3][eigensystem[3]];
    inner=abs(innerProduct(eigenvector, states, eigenvector2, states2));
    return inner;
end

function calculateFidelityPlain(hmin, hmax, num, N, J, J2, bonds)
    println("begin");
    println(J);
    fidelities=Any[];
    hstep=hmax/num-hmin/num;
    all=Any[];
    println("hlist begin");
    hlist=generateHListUniform(hmin, hmax, num);
    println("fidelity list begin");
    list=generateFidelityList(hmin, hmax, num, N, J, J2, bonds);
    println("fidelity list length", length(list));
    eigenvectors=list[1];
    states=list[2];
    for i=1:length(eigenvectors)-1
        inner=abs(innerProduct(eigenvectors[i], states[i], eigenvectors[i+1], states[i+1]));
        push!(fidelities,inner);
    end
    return fidelities;
end

function generateFidelityListNoSymmetry(hmin, hmax, num, N, J, J2, bonds)
    hstep=hmax/num-hmin/num;
    i=hmin;
    eigenvectors=Any[];
    states=Any[];
    all=Any[];
    println(J);

    while(i<=hmax+hstep)
        println("on i: ", i);
        temp=calculateEigensystemTransverseNoSymmetry(
                num,
                J,
                J2,
                hs,
                hs2,
                bonds,
                eigmethod,
                num,
                "H1",
            )
            temp=calculateEigensystemTransverseNoSymmetry(N, J, J2, i, bonds,"lanczos", "one", i, 0, "H1");

        #eigensystem=getLowestLyingStates(temp[1], temp[2]);
        eigensystem=getLowestLyingStates(temp[1], temp[2]);
        eigenvector=eigensystem[2];
        state=temp[3][eigensystem[3]];
        push!(eigenvectors, eigenvector);
        push!(states, state);
        if(hstep<=0)
            break;
        end
        i+=hstep;
    end
    push!(all, eigenvectors);
    push!(all, states);
    return all;
end

function calculate_correlation_length(eigenvector,r,n)
    sum=0
    for i=0:((2^n)-1)
        spin1=getTi(0,i) == 1 ? 1 : -1
        spin2=getTi(r-1,i) ==1 ? 1 : -1
        sum+=spin1*spin2*abs2(eigenvector[i+1])
    end
    return sum
end

function calculateSzSquared(eigenvector, N)
    sum=0;
    for j=1:N
        for k=1:N
            for i=0:2^(N)-1
                spin1=getTi(j-1,i) == 1 ? 1 : -1
                spin2=getTi(k-1,i) == 1 ? 1 : -1
                sum+=spin1*spin2*abs2(eigenvector[i+1])
            end
    end
end
    return sum/N^2
end


function calculate_participation_ratio(eigenvector)
    sum1=0
    sum2=0
    for i=1:length(eigenvector)
        temp=eigenvector[i]
        sum1+=abs(temp)^4
        sum2+=abs(temp)^2
    end
    return sum1/(sum2)^2
end

function expand_energy_eigenbasis(eigenvectors, state, states)
    coefficients=[]
    indicies=[]
    #println(state)
    for i=1:size(eigenvectors,2)
        in=innerProduct(state, eigenvectors[:,i],states)
        push!(coefficients,in)

        if(abs(in)>10^(-5))
            push!(indicies, i)
        end
    end
    return (indicies,coefficients)
end

function calculate_Time_Evolved_X(state, eigenvalues, eigenvectors, N, t, states, o)
    sum=0;
    for i=1:length(state)
        sum+=abs2(state[i])*get_average_ith_x(eigenvectors[:,i], states, N, o)
        for j=i+1:length(state)
            if(j!=i)
                sum+=conj(state[i])*state[j]*exp(-im*(eigenvalues[i]-eigenvalues[j])*t)*get_average_ith_x(eigenvectors[:,i],eigenvectors[:,j], states, N,o)+conj(state[j])*state[i]*exp(im*(eigenvalues[i]-eigenvalues[j])*t)*get_average_ith_x(eigenvectors[:,i],eigenvectors[:,j], states, N,o)
            end
        end
    end
    return sum
end


function calculate_Time_Evolved_Z_Pre(state, eigenvectors, N, o)
    sum=0;
    for i=1:length(state)
        sum+=abs2(state[i])*get_average_ith_z(eigenvectors[:,i], N, o)
    end
    return sum
end

function calculate_Time_Evolved_Z_Pre(state, eigenvectors, N, o, matrix)
    sum=0;
    for i=1:length(state)
        sum+=abs2(state[i])*matrix[i,i]
    end
    return sum
end

function calculate_Time_Evolved_Z(state, eigenvalues, eigenvectors, N, t, states, o, matrix)
    sum=0;
    for i=1:length(state)
        for j=i+1:length(state)
            #sum+=2*state[i]*state[j]*cos((eigenvalues[i]-eigenvalues[j])*t)*(matrix[j,i]+matrix[i,j])
            one=conj(state[i])*state[j]*exp(im*(eigenvalues[j]-eigenvalues[i])*t)*matrix[i,j]
            two=conj(state[j])*state[i]*exp(im*(eigenvalues[i]-eigenvalues[j])*t)*matrix[j,i]
            sum+=one+two
            #println(sum)
            #sum+=(conj(state[i])*state[j]*exp(im*(eigenvalues[i]-eigenvalues[j])*t)*matrix[j,i]+conj(state[j])*state[i]*exp(-im*(eigenvalues[i]-eigenvalues[j])*t)*matrix[i, j])
        end
    end
    return sum
end

function calculate_Time_Evolved_Y(state, eigenvalues, eigenvectors, N, t)
    con=0;
    for i=1:length(state)
        con+=abs2(state[i])*calculateSy(eigenvectors[i],eigenvectors[i], N)
    end
    sum=0
    for i=1:length(state)
        for j=1:length(state)
            if(j!=i)
                sum+=conj(state[i])*state[m]*exp(-im*eigenvalues[i]-eigenvalues[j]*t)*calculateSy(eigenvectors[i],eigenvectors[j], N)
            end
        end
    end
    return sum+con
end

# state: list of all the complex coefficients in front of each eigenvector
function time_evolve_state(state, eigenvalues, t)
    state_temp = copy(state)
    for i=1:length(state)
        state_temp[i] = exp(-im*eigenvalues[i]*t) * state[i]
    end
    return state_temp
end

# use the same initial state as in the PRX paper
function generate_initial_state(theta, N)
    a = cos(theta / 2)
    b = im * sin(theta / 2)
    psi0 = [a, b]

    state = psi0
    for i=2:N        
        state = kron(state, psi0)
    end
    return state
end

# get the probability of measuring the nth computational basis state, given a state vector (2^N list of coefficients in the
# superposition of eigenvectors) and all the eigenstates
# TASK: calculate the probability of measuring a random computational basis state n, over time, over many disorder realizations

function get_probability_of_n_computational_basis(n, states, eigenstates)
    N=size(eigenstates, 2)
    temp=0
    for i=1:length(states)
        temp+=states[i]*eigenstates[n, i]
    end
    return abs2(temp)
end

function calculateSz(eigenvector, eigenvector2, N)
    sum=0;
    for i=0:2^(N)-1
        temp=countBits(i);
        sum+=(temp-(N-temp))*conj(eigenvector2[i+1])*eigenvector[i+1];
    end
    return (1/(N))*(1/2)*sum;
end



function calculateSy(eigenvector, eigenvector2, N)
    sum=0;
    for i=0:2^(N)-1
        for j=1:N
            b=flipBit(j-1, i);
            factor=getTi(j-1, i)==1 ? im : -im
            sum+=conj(eigenvector2[b+1])*factor*eigenvector[i+1]
        end
    end
    return sum;
end


function calculate_thermal_state_energy(eigenvalues, beta)
    sum=0
    part=0
    for i=1:length(eigenvalues)
        sum+=exp(-beta*eigenvalues[i])*eigenvalues[i]
        part+=exp(-beta*eigenvalues[i])
    end
    return sum/part
end


function calculate_thermal_state_x(eigenvalues, eigenvectors, beta, N)
    sum=0
    part=0
    for i=1:length(eigenvalues)
        sum+=exp(-beta*eigenvalues[i])*calculateSx(eigenvectors[:,i], N)
        part+=exp(-beta*eigenvalues[i])
    end
    return sum/part
end

function calculate_thermal_state_z(eigenvalues, eigenvectors, beta, N)
    sum=0
    part=0
    for i=1:length(eigenvalues)
        sum+=exp(-beta*eigenvalues[i])*calculateSz(eigenvectors[:,i], N)
        part+=exp(-beta*eigenvalues[i])
    end
    return sum/part
end

function calculate_thermal_state_y(eigenvalues, eigenvectors, beta, N)
    sum=0
    part=0
    for i=1:length(eigenvalues)
        sum+=exp(-beta*eigenvalues[i])*calculateSy(eigenvectors[:,i],eigenvectors[:,i], N)
        part+=exp(-beta*eigenvalues[i])
    end
    return sum/part
end

function isMomentumZero(state1, state2, N)
    states=[]
    states2=[]
    for i=0:2^N-1
        push!(states, i)
    end
    for i=1:length(state1)
        temp=circular_shift(states[i], N,1)
        push!(states2,temp)
    end
    #println(abs(innerProduct(state1, states,state2,states2)))
    return (abs(innerProduct(state1, states,state2,states2)-1))<0.000001
end

function get_average_ith_z(eigenvector,N, j)
    sum=0
    for i=1:length(eigenvector)
        spin=getTi(j, eigenvector[i])==1 ? 1 : -1
        sum+=abs2(eigenvector[i])*spin
    end
    return sum
end

function find_which_eigenvectors(state,eigenvectors, states)
    temp=[]
    for i in range(1,length(eigenvectors[:,1]))
        if(innerProduct(state, eigenvectors[:,i],states)>(1/length(states)))
            push!(temp, i)
        end
    end
    return temp
end

function calculate_matrix_elements_j_z(eigenvectors, N, j)
    temp=zeros(Complex, size(eigenvectors,2),  size(eigenvectors,2))
    Threads.@threads for i in range(1, size(eigenvectors,2))
        for k in range(1,  size(eigenvectors,2))
            temp[i, k]=get_average_ith_z(eigenvectors[:,i], eigenvectors[:,k], N, j)
        end
    end
    return temp
end


function calculate_matrix_elements_entanglement(eigenvectors, N, states, listA)
    temp=zeros(Complex, size(eigenvectors,2),  size(eigenvectors,2))
    for i in range(1, size(eigenvectors,2))
        for k in range(1,  size(eigenvectors,2))
            temp[i, k]= getEntanglementEntropy(eigenvectors[:,i], states,listA, N)
        end
    end
    return temp
end


function calculate_matrix_elements_j_x(eigenvectors, N, j)
    temp=zeros(Complex, size(eigenvectors,2),  size(eigenvectors,2))
    for i in range(1, size(eigenvectors,2))
        for k in range(1, size(eigenvectors,2))
            temp[i, k]=get_average_ith_x(eigenvectors[:,i], eigenvectors[:,k], N, j)
            #println("i="*string(i))
            #println("k="*string(k))
        end
    end
    #println(ishermitian(temp))
    return temp
end

function calculate_matrix_elements_j_y(eigenvectors, N, j)
    temp=zeros(Complex, size(eigenvectors,2),  size(eigenvectors,2))
    for i in range(1, size(eigenvectors,2))
        for k in range(1, size(eigenvectors,2))
            temp[i, k]=get_average_ith_y(eigenvectors[:,i], eigenvectors[:,k], N, j)
            #println("i="*string(i))
            #println("k="*string(k))
        end
    end
    #println(ishermitian(temp))
    return temp
end

function get_average_ith_z(eigenvector, eigenvector2,N, j)
    sum=0
    for i=1:length(eigenvector[:,1])
        spin=getTi(j-1, i-1)==1 ? 1 : -1
        sum+=conj(eigenvector2[i])*(eigenvector[i])*spin
    end
    return sum
end


function get_average_ith_x(eigenvector, N, j)
    sum=0;
    for i=0:2^(N)-1
        b=flipBit(j-1, i);
        sum+=eigenvector[i+1]*conj(eigenvector[b+1]);
    end
    return sum;
end


function get_average_ith_x(eigenvector, eigenvector2, N, j)
    sum=0;
    for i=0:2^(N)-1
        b=flipBit(j-1, i);
        sum+=eigenvector[i+1]*conj(eigenvector2[b+1]);
    end
    return sum;
end


function get_average_ith_y(eigenvector, N, j)
    sum=0;
    for i=0:2^(N)-1
        b=flipBit(j-1, i);
        factor=getTi(j-1, i)==1 ?  im : -im
        sum+=factor*eigenvector[i+1]*conj(eigenvector[b+1]);
    end
    return sum;
end



function get_average_ith_y(eigenvector, eigenvector2, N, j)
    sum=0;
    for i=0:2^(N)-1
        b=flipBit(j-1, i);
        factor=getTi(j-1, i)==1 ?  im : -im
        sum+=factor*eigenvector[i+1]*conj(eigenvector2[b+1]);
    end
    return sum;
end

function get_time_evolved_state(state, eigenvalues, eigenvectors, t)
    temp = zeros(Float64, length(eigenvectors[:, 1]))
    for i=1:length(eigenvalues)
        temp += (state[i] * exp(-im * t * eigenvalues[i])) .* eigenvectors[:, i]
    end
    return temp
end