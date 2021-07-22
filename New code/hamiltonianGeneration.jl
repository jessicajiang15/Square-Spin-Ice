include("bondGeneration.jl");
include("momentumHelpers.jl")
include("reflectionHelper.jl")
using SparseArrays
function constructTransverseHamiltonian(states, bonds, N, J, eigmethod, randList)
    #loop through ALL the possible states...
    list=states[1];
    map=states[2];
    #println(states);
    local H;
    if(eigmethod=="full")
        H=zeros(Float64, length(list),length(list));
    else
        H=spzeros(Float64, length(list),length(list));
    end
    #H::Matrix{Float64}=zeros(Int, length(list),length(list));

    for i=1:length(list)
        #NOW loop through all the possible SITES
        count=0;
        for j=1:N*N
            #for this particular site, find ALL of its bonds
            #println("bonds length", length(bonds));
            theThing= getTi(j,list[i]) == 1 ? randList[j]/2 : -randList[j]/2;
            #1/2 for double counting
            #println((1/2)*theThing);
            H[i,i]+= theThing;
            for z=1:length(bonds)
                #get the number of the thing it is bonding with!
                if containsSite(j,bonds[z])
                    bond1::Int=bonds[z].site1.num;
                    bond2::Int=bonds[z].site2.num;
                    if(bond1!=j)
                        temp::Int=bond1;
                        bond1=bond2;
                        bond2=temp;
                    end
                        #flip the bits at those relevant places
                        b::Int=flipBits(bond1-1, bond2-1,list[i]);
                        t::Int=map[b];
                        H[i,t]=J/4;
                        end
                    end

                end
                #println("i count ", i, count);
            end
        return H;
end

function constructTransverseHamiltonianNoSymmetrySx(bonds, N, J, eigmethod, randomList)
    if(eigmethod=="full")
        H=zeros(Float64, 2^(N*N), 2^(N*N));
    else
        H=spzeros(Float64, 2^(N*N), 2^(N*N));
    end
    #H::Matrix{Float64}=zeros(Int, length(list),length(list));
    for i=0:2^(N*N)-1
        #NOW loop through all the possible SITES
        for j=1:N*N
            b::Int=flipBit(j, i);
            H[i+1,b+1]=randList[j]/2;
            #for this particular site, find ALL of its bonds
            for z=1:length(bonds)
                #get the number of the thing it is bonding with!
                if containsSite(j,bonds[z])
                    bond1::Int=bonds[z].site1.num;
                    bond2::Int=bonds[z].site2.num;
                    local theThing;
                    if(getTi(bond1-1,i)==getTi(bond2-1,i))
                        theThing=J/4;
                    else
                        theThing=-J/4;
                    end
                    H[i+1,i+1]+= (1/2)*theThing;
                        end
                    end
                end
            end
        return H;
end


function constructTransverseHamiltonianNoSymmetrySz(states, bonds, N, J, eigmethod, randList)
    #loop through ALL the possible states...
    #println(states);
    local H;
    if(eigmethod=="full")
        H=zeros(Float64, 2^(N*N), 2^(N*N));
    else
        H=spzeros(Float64,  2^(N*N), 2^(N*N));
    end
    #H::Matrix{Float64}=zeros(Int, length(list),length(list));

    for i=0:2^(N*N)-1
        #NOW loop through all the possible SITES
        count=0;
        for j=1:N*N
            #for this particular site, find ALL of its bonds
            #println("bonds length", length(bonds));
            theThing= getTi(j,i) == 1 ? randList[j]/2 : -randList[j]/2;
            #1/2 for double counting
            #println((1/2)*theThing);
            H[i+1,i+1]+= theThing;
            for z=1:length(bonds)
                #get the number of the thing it is bonding with!
                if containsSite(j,bonds[z])
                    bond1::Int=bonds[z].site1.num;
                    bond2::Int=bonds[z].site2.num;
                    if(bond1!=j)
                        temp::Int=bond1;
                        bond1=bond2;
                        bond2=temp;
                    end
                        #flip the bits at those relevant places
                        b::Int=flipBits(bond1-1, bond2-1,i);
                        H[i+1,b+1]=J/4;
                        end
                    end

                end
                #println("i count ", i, count);
            end
        return H;
end

#TODO: update this and update the other one
function constructTransverseHamiltonianSzBasis(states, bonds, N, J, eigmethod, randList)
    #loop through ALL the possible states...
    list=states[1];
    map=states[2];
    #println(states);
    local H;
    if(eigmethod=="full")
        H=zeros(Float64, length(list),length(list));
    else
        H=spzeros(Float64, length(list),length(list));
    end
    #H::Matrix{Float64}=zeros(Int, length(list),length(list));
    for i=1:length(list)
        #NOW loop through all the possible SITES
        for j=1:N*N
            #for this particular site, find ALL of its bonds
            for z=1:length(bonds)
                #get the number of the thing it is bonding with!
                if containsSite(j,bonds[z])
                    bond1::Int=bonds[z].site1.num;
                    bond2::Int=bonds[z].site2.num;
                    local theThing;
                    if(getTi(bond1-1,list[i])==getTi(bond2-1,list[i]))
                        theThing=randList[z]/4;
                    else
                        theThing=-randList[z]/4;
                    end
                    H[i,i]+= (1/2)*theThing;
                    #1/2 for double counting
                        #flip the bits at those relevant places
                        b::Int=flipBit(bond1-1, list[i]);
                        t::Int=map[b];
                        H[i,t]=J/4;
                        end
                    end
                end
            end
        return H;
end


function constructHeisenbergHamiltonian(states, bonds, N, J, eigmethod)
    #loop through ALL the possible states...
    list=states[1];
    map=states[2];
    #println(states);
    local H;
    if(eigmethod=="full")
        H=zeros(Float64, length(list),length(list));
    else
        H=spzeros(Float64, length(list),length(list));
        println(typeof(H));
    end    #H::Matrix{Float64}=zeros(Int, length(list),length(list));
    for i=1:length(list)
        #NOW loop through all the possible SITES
        for j=1:N*N
            #for this particular site, find ALL of its bonds
            for z=1:length(bonds)
                #get the number of the thing it is bonding with!
                if containsSite(j,bonds[z])
                    bond1::Int=bonds[z].site1.num;
                    bond2::Int=bonds[z].site2.num;
                    #now, are the spins in those two places the same? if so, the
                    #diagonal entry at i,i is 1
                    if(getTi(bond1-1,list[i])==getTi(bond2-1, list[i]))
                        H[i, i]+=(1/2)*J/4;
                    else
                        H[i, i]-=(1/2)*J/4;
                        #flip the bits at those relevant places
                        b::Int=flipBits(bond1-1, bond2-1,list[i]);
                        t::Int=map[b];
                            H[i,t]=J/2;
                    end

                end
            end
        end
    end
    states=0;
    GC.gc();
            return H;
end

#refstates map tells u which refstate is the same
function constructHamiltonianHeisenbergMomentum(refStates, N, psector, bonds, refStatesMap, eigmethod)
#array of viable ref states
    ref::Array{refState}=refStates[1];
    #map that maps xstates to their numbering within hamiltonian
    map::Dict{Int, Int}=refStates[2];
    local H;
    if(eigmethod=="full")
        H=zeros(Complex{Float64}, length(ref),length(ref));
    else
        H=spzeros(Complex{Float64}, length(ref),length(ref));
    end
    println(size(H));
    for i=1:length(ref)
        #loop over all sites
            for j=1:N*N
                #check which bonds contain that site
                for z=1:length(bonds)
                    #if it does we proceed
                    if(containsSite(j, bonds[z]))
                        bond1::Int=bonds[z].site1.num;
                        bond2::Int=bonds[z].site2.num;
                        if(bond1!=j)
                            temp::Int=bond1;
                            bond1=bond2;
                            bond2=temp;
                        end
                        if(getTi(bond1-1, ref[i].state)==getTi(bond2-1, ref[i].state))
                            H[i,i]+=(1/2)*1/4;
                        else
                            H[i,i]-=(1/2)*1/4;
                            b::Int=flipBits(bond1-1, bond2-1, ref[i].state);
                            #c::Int=-1;
                            #p, N, ref
                            if(isViable(psector, N, refStatesMap[b].ref))
                                c=map[refStatesMap[b].ref.state];
                                norm::Float64=sqrt(ref[i].periodicity/refStatesMap[b].ref.periodicity);
                                #norm::Float64=sqrt((N^2)/ref[i].periodicity))/(((((N^2)/refStatesMap[b].ref.periodicity))));
                                H[i, c]+=(1/2)*exp(-1im*calculateMomentum(psector, N)*refStatesMap[b].shiftsNeeded)*(1/2)*norm;
                            end
                        end
                    end
                end
            end
    end
    return H;

end



#refstates map tells u which refstate is the same
function constructHamiltonianHeisenbergMomentum2d(refStates, N, momentum, bonds, refStatesMap, eigmethod)
#array of viable ref states
    ref::Array{refState2d}=refStates[1];
    phaseFactorSums::Dict{Int, Complex{Float64}}=refStates[3];
    #map that maps xstates to their numbering within hamiltonian
    map::Dict{Int, Int}=refStates[2];
    local H;
    if(eigmethod=="full")
        H=zeros(Complex{Float64}, length(ref),length(ref));
    else
        H=spzeros(Complex{Float64}, length(ref),length(ref));
    end
        println(size(H));
    for i=1:length(ref)
        #loop over all sites
            for j=1:N*N
                #check which bonds contain that site
                for z=1:length(bonds)
                    #if it does we proceed
                    if(containsSite(j, bonds[z]))
                        bond1::Int=bonds[z].site1.num;
                        bond2::Int=bonds[z].site2.num;
                        if(getTi(bond1-1, ref[i].state)==getTi(bond2-1, ref[i].state))
                            H[i,i]+=(1/2)*1/4;
                        else
                            H[i,i]-=(1/2)*1/4;
                            b::Int=flipBits(bond1-1, bond2-1, ref[i].state);
                            #c::Int=-1;
                            #p, N, ref
                            if(refStatesMap[b].ref.state in keys(map))
                                local c=map[refStatesMap[b].ref.state];
                                one=ref[i].numUniqueSt;
                                two=phaseFactorSums[ref[i].state];
                                three=refStatesMap[b].ref.numUniqueSt;
                                four=phaseFactorSums[refStatesMap[b].ref.state];
                                Na=one*abs2(two);
                                Nb=three*abs2(four);
                                #println("Na, ", Na, ", Nb: ", Nb);
                                norm::Float64=sqrt(Nb/Na);
                                #work on this
                                H[i, c]+=(1/2)*exp(-1im*(calculateMomentum(momentum.px, N)*refStatesMap[b].shiftsXNeeded+calculateMomentum(momentum.py, N)*refStatesMap[b].shiftsYNeeded))*(1/2)*norm;
                            end
                        end
                    end
                end
            end
    end
    println("h00 ", H[1,1]);
    return H;

end




#refstates map tells u which refstate is the same
function constructHamiltonianHeisenbergReflection(refStates, N, reflection, bonds, refStatesMap, eigmethod)
#array of viable ref states
    ref::Array{reflectionRefState}=refStates[1];
    #map that maps xstates to their numbering within hamiltonian
    map::Dict{Int, Int}=refStates[2];
    local H;
    if(eigmethod=="full")
        H=zeros(Float64, length(ref),length(ref));
    else
        H=spzeros(Float64, length(ref),length(ref));
    end
        println(size(H));
    for i=1:length(ref)
        #loop over all sites
            for j=1:N*N
                #check which bonds contain that site
                for z=1:length(bonds)
                    #if it does we proceed
                    if(containsSite(j, bonds[z]))
                        bond1::Int=bonds[z].site1.num;
                        bond2::Int=bonds[z].site2.num;
                        if(bond1!=j)
                            temp::Int=bond1;
                            bond1=bond2;
                            bond2=temp;
                        end
                        if(getTi(bond1-1, ref[i].state)==getTi(bond2-1, ref[i].state))
                            H[i,i]+=(1/2)*1/4;
                        else
                            H[i,i]-=(1/2)*1/4;
                            b::Int=flipBits(bond1-1, bond2-1, ref[i].state);
                            #c::Int=-1;
                            #p, N, ref
                            if(refStatesMap[b].refState.state in keys(map))
                                local c=map[refStatesMap[b].refState.state];
                                norm::Float64=sqrt(refStatesMap[b].refState.uniqueReflections/ref[i].uniqueReflections);
                                H[i, c]+=(1/2)*(refStatesMap[b].refXNeeded == 1 ? refStatesMap[b].refState.xC : 1)*(refStatesMap[b].refYNeeded == 1 ? refStatesMap[b].refState.yC : 1)*norm*(1/2);
                                #H[i, c]+=(1/2)*(refStatesMap[b].refState.xC)*(refStatesMap[b].refState.yC)*norm*(1/2);
                            end
                        end
                    end
                end
            end
    end
    return H;

end
