include("bondGeneration.jl");
include("momentumHelpers.jl")
using SparseArrays
function constructTransverseHamiltonian(states, bonds, N, J, h)
    #loop through ALL the possible states...
    list=states[1];
    map=states[2];
    #println(states);
    H::SparseMatrixCSC{Float64}=spzeros(Int, length(list),length(list));
    for i=1:length(list)
        #NOW loop through all the possible SITES
        for j=1:N*N
            #for this particular site, find ALL of its bonds
            for z=1:length(bonds)
                #get the number of the thing it is bonding with!
                if containsSite(j,bonds[z])
                    bond1::Int=bonds[z].site1.num;
                    bond2::Int=bonds[z].site2.num;
                    theThing= getTi(bond1-1,list[i]) == 1 ? h/2 : -h/2;
                    H[i,i]+= theThing;
                        #flip the bits at those relevant places
                        b::Int=flipBits(bond1-1, bond2-1,list[i]);
                        t::Int=map[b];
                        H[i,t]=J/4;
                        end
                    end

                end
            end
        return H;
end


function constructHeisenbergHamiltonian(states, bonds, N, J)
    #loop through ALL the possible states...
    list=states[1];
    map=states[2];
    #println(states);
    H::SparseMatrixCSC{Float64}=spzeros(Int, length(list),length(list));
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
                    if(bond1!=j)
                        temp::Int=bond1;
                        bond1=bond2;
                        bond2=temp;
                    end
                    if(getTi(bond1-1,list[i])==getTi(bond2-1, list[i]))
                        H[i, i]+=(1/2)*J/4;
                    else
                        H[i, i]-=(1/2)*J/4;
                        #flip the bits at those relevant places
                        b::Int=flipBits(bond1-1, bond2-1,list[i]);
                        t::Int=map[b];
                            H[i,t]=(1/2)*J/2;
                    end

                end
            end
        end
    end
            return H;
end

#refstates map tells u which refstate is the same
function constructHamiltonianHeisenbergMomentum(refStates, N, psector, bonds, refStatesMap)
#array of viable ref states
    ref::Array{refState}=refStates[1];
    #map that maps xstates to their numbering within hamiltonian
    map::Dict{Int, Int}=refStates[2];
    H::SparseMatrixCSC{Float64}=spzeros(Float64, length(ref), length(ref));
    println(size(H));
    for i=1:length(ref)
            for j=1:N*N
                for z=1:length(bonds)
                    if(containsSite(j, bonds[z]))
                        bond1::Int=bonds[z].site1.num;
                        bond2::Int=bonds[z].site2.num;
                        if(bond1!=j)
                            temp::Int=bond1;
                            bond1=bond2;
                            bond2=temp;
                        end
                        if(getTi(bond1-1, ref[i].state)==getTi(bond2-1, ref[i].state))
                            H[i,i]+=1/4;
                        else
                            H[i,i]-=1/4;
                        end
                        b::Int=flipBits(bond1-1, bond2-1, ref[i].state);
                        c::Int=-1;
                        try
                        c=map[refStatesMap[b].state];
                        catch
                        continue;
                    end
                        norm::Float64=sqrt(ref[i].periodicity/ref[c].periodicity);
                        H[i, c]+=exp(1im*2*pi*calculateMomentum(psector, N)*refStatesMap[b].shiftsNeeded/(N*N))+(1/2)*norm;
                    end
                end
            end
    end
    return H;

end
