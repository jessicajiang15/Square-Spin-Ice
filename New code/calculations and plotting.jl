using Plots
include("necessaryStructures.jl")

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


function calculateSz(eigenvector, N)
    sum=0;
    for i=0:2^(N*N)-1
        temp=countBits(i);
        sum+=(temp-(N*N-temp))*abs2(eigenvector[i+1]);
    end
    return (1/(N*N))*(1/2)*sum;
end

function calculateSz(eigenvector, states, N)
    sum=0;
    for i=1:length(states)
        temp=countBits(states[i]);
        sum+=(1/2)*(temp-(N*N-temp))*conj(eigenvector[i])*eigenvector[i];
    end
    return -(1/(N*N))*sum;
end

function szMatrix(eigenvector, states, N)
    m::Matrix{Float64}=zeros(length(eigenvector), length(eigenvector));
    for i=length(states)
        temp=countBits(states[i]);
        m[i, i]+=(1/2)*(temp-(N*N-temp));
    end
    return conj.(eigenvector')*m*eigenvector*(1/(N*N));
end

function calculateSx(eigenvector, N)
    sum=0;
    for i=0:2^(N*N)-1
        for i=1:N*N
            b=flipBit(i-1, i);
            sum+=eigenvector[i+1]*conj(eigenvector[b+1]);
        end
    end
    return (1/2)*sum;
end

#let state be map from state -> position in array
function calculateSx(eigenvector, states, map, N)
    sum=0;
    for i=1:length(states)
        for j=1:N*N
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
    one=getTi(state, starting);
    two=getTi(state, starting+1);
    three=getTi(state, starting+N);
    four=getTi(state, starting+1+N);
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
    println("length eigenvalues, ", length(eigenvalues));
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
        index=binarySearchIt(states[i], states2);
        #println("states[i]", states[i]);
        #println("index", index);
        if(index!=-1)
            sum+=conj(eigenvector[i])*eigenvector2[index];
        end
    end
    return sum;
end

#assumes the two states are the same
function innerProduct(eigenvector, states, eigenvector2)
    sum=0;
    for i=1:length(states)
           sum+=conj(eigenvector[i])*eigenvector2[i];
    end
    return sum;
end

#please just assume that the list of states will be the same for the new eigenvector and stuff...
function calculateFidelity(eigenvector, states, h, N, deltah, bonds, J)
    temp=calculateEigensystemTransverse(N, J, h+deltah, bonds,"lanczos", "one", h+deltah, 0);
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


function calculateFidelity(eigenvector, states, h, N, deltah, bonds, J)
    temp=calculateEigensystemTransverse(N, J, h+deltah, bonds,"lanczos", "one", h+deltah, 0);
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

function generateFidelityList(hmin, hmax, num, N, J, bonds)
    hstep=hmax/num-hmin/num;
    i=hmin;
    eigenvectors=Any[];
    states=Any[];
    all=Any[];
    println(J);

    while(i<=hmax+hstep)
        println("on i: ", i);
        temp=calculateEigensystemTransverse(N, J, i, bonds,"lanczos", "one", i, 0);
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

function calculateFidelity(hmin, hmax, num, N, J, bonds)
    println("begin");
    println(J);
    fidelities=Any[];
    hstep=hmax/num-hmin/num;
    all=Any[];
    println("hlist begin");
    hlist=generateHListUniform(hmin, hmax, num);
    println("fidelity list begin");
    list=generateFidelityList(hmin, hmax, num, N, J, bonds);
    println("fidelity list length", length(list));
    eigenvectors=list[1];
    states=list[2];
    for i=1:length(eigenvectors)-1
        inner=abs(innerProduct(eigenvectors[i], states[i], eigenvectors[i+1], states[i+1]));
        push!(fidelities,2*(1-inner)/(hstep^2));
    end
    return fidelities;
end

function getSzPlaquette(state, plaquetteIndicies)
    sum=0;
    for i=1:length(plaquetteIndicies)
        sum+=getTi(plaquetteIndicies(i), state)
    end
    return sum;
end

function calculateSxIndicies(eigenvector, states, map, N, indicies)
    sum=0;
    for i=1:length(states)
        for j=1:length(indicies)
            b=flipBit(indicies[j]-1, states[i]);
            if(b in keys(map))
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
                jSz=getSzPlaquette(state, plaquettes[j]);
                zSz=getSzPlaquette(state, plaquettes[z]);
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
