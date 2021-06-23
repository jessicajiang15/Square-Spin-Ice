
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


function calcualteExpESq(beta, E)
    sum=0;
    for i=1:length(E)
        sum+=(E[i]^2)*exp(E[i]*beta);
    end
    sum/=partition(beta, E);
    return sum;
end

function calculateHeatCapacity(beta, E, E2)
    return (beta^2)*(E2-E^2);
end

function calculateG(beta, partition)
    return (-1/beta)*log(partition);
end

function calculateEntropy(beta, G, E2)
    return beta*(G-H)
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
    heatCapcities=Float64[];
    gs=Float64[];
    expESqs=Float64[];
    expEs=Float64[];
    partitions=Float64[];
    entropies=Float64[];
    betas=Float64[];
    b=bmin;

    while b<=bmax
        partition=partition(b, eigenvalues, partition);
        push!(partitions, partition);
        expE=calculateExpE(b, eigenvalues, partition);
        push!(expEs, expE);
        expESq=calculateExpESq(b, eigenvalues, partition);
        push!(expESqs, expESq);
        heatCap=calculateHeatCapacity(b, expE, expESq);
        push!(heatCapacities, heatCap);
        g=calculateG(b, partition);
        push!(heatCapacities, heatCap);
        entropy=calculateEntropy(beta, g, expESq);
        push!(entropies, entropy);
        b+=bstep;
    end
    push!(all, betas);
    push!(all, partitions);
    push!(all, expEs);
    push!(all, expESqs);
    push!(all, heatCapcities);
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
function theplot(list, i)
    plot(list[1], list[i]);
end

function plotQuantity(bmin, bmax, bstep, eigenvalues, i)
    all=getAllRelevantQuantities(bmin, bmax, bstep, eigenvalues);
    str::String;
    if(i==5)
        str="Heat Capacity";
    elseif(i==6)
        str="Free Energy";
    else
        str="Entropy";
    end
    theplot(all, i, title = string(str, "vs. 1/kbT"), label = [str]);
end


function plotQuantity(all, i)
    str::String;
    if(i==5)
        str="Heat Capacity";
    elseif(i==6)
        str="Free Energy";
    else
        str="Entropy";
    end
    theplot(all, i, title = string(str, "vs. 1/kbT"), label = [str]);
end
