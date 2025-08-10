include("bondGeneration.jl");
include("momentumHelpers.jl")
include("reflectionHelper.jl")
include("andersonHelper.jl")
include("renormalizationGroupHelpers.jl")
using SparseArrays
#IMPORTANT: N IS THE TOTAL NUMBER OF SITES!!!!!!!!!!!!
function constructTransverseHamiltonian(states, bonds, N, J, J2, eigmethod, randList)
    #loop through ALL the possible states...
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];

    list=states[1];
    map=states[2];
    #println(states);
    local H;
    if(eigmethod=="full")
        H=zeros(Float64, length(list),length(list));
        #=
    else
    H=spzeros(Float64, length(list),length(list));
    =#
end
#H::Matrix{Float64}=zeros(Int, length(list),length(list));

for i=1:length(list)
    #NOW loop through all the possible SITES
    count=0;
    for j=1:N
        #for this particular site, find ALL of its bonds
        #println("bonds length", length(bonds));
        theThing= getTi(j-1,list[i]) == 1 ? randList[j]/2 : -randList[j]/2;
        #println((1/2)*theThing);
        push!(rows, i);
        push!(cols, i);
        push!(values, theThing);
        if(eigmethod=="full")
            H[i,i]+=theThing;
        end
        #
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
                a= bonds[z].isNear ? J : J2;
                if(eigmethod=="full")
                    H[i,t]=a/4;
                end
                push!(rows, i);
                push!(cols, t);
                push!(values, (1/2)*a/4);
            end
        end

    end
    #println("i count ", i, count);
end
if(eigmethod!="full")
    H=sparse(rows, cols, values);
end
return H;
end

#H1
function constructTransverseHamiltonianNoSymmetrySx(bonds, N, J, J2, eigmethod, randList)
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];

    if(eigmethod=="full")
        H=zeros(Float64, 2^(N), 2^(N));
    end
    #H::Matrix{Float64}=zeros(Int, length(list),length(list));
    for i=0:2^(N)-1
        #NOW loop through all the possible SITES
        for j=1:N
            b::Int=flipBit(j-1, i);
            push!(rows, i+1);
            push!(cols, b+1);
            push!(values, randList[j]/2);
            if(eigmethod=="full")
                H[i+1,b+1]=randList[j]/2;
            end
            #for this particular site, find ALL of its bonds
            for z=1:length(bonds)
                #get the number of the thing it is bonding with!
                if containsSite(j,bonds[z])
                    bond1::Int=bonds[z].site1.num;
                    bond2::Int=bonds[z].site2.num;
                    local theThing;
                    a= bonds[z].isNear ? J : J2;
                    #println("bond1", bond1);
                    #println("bond2", bond2);

                    if(getTi(bond1-1,i)==getTi(bond2-1,i))
                        theThing=a/4;
                    else
                        theThing=-a/4;
                    end
                    push!(rows, i+1);
                    push!(cols, i+1);
                    push!(values, (1/2)*theThing);
                    if(eigmethod=="full")
                        H[i+1,i+1]+= (1/2)*theThing;
                    end
                end
            end
        end
    end
    if(eigmethod!="full")
        H=sparse(rows, cols, values);
    end
    if(N<=1)
        println("num sites, ", N, ", H", Matrix(H))
    end
    return H;
end

#H2, this is sxsx+sz so field is in the z direction
function constructTransverseHamiltonianNoSymmetrySz(bonds, N, J, J2, eigmethod, randList)
    #loop through ALL the possible states...
    #println(states);
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];
    local H;

    if(eigmethod=="full")
        H=zeros(Float64, 2^(N), 2^(N));
    end
    #H::Matrix{Float64}=zeros(Int, length(list),length(list));

    for i=0:2^(N)-1
        #NOW loop through all the possible SITES
        count=0;
        for j=1:N
            #for this particular site, find ALL of its bonds
            #println("bonds length", length(bonds));
            theThing= getTi(j-1,i) == 1 ? randList[j]/2 : -randList[j]/2;
            #1/2 for double counting
            #println((1/2)*theThing);
            push!(cols, i+1);
            push!(rows, i+1);
            push!(values, theThing);
            if(eigmethod=="full")
                H[i+1,i+1]+= theThing;
            end
            for z=1:length(bonds)
                #get the number of the thing it is bonding with!
                if containsSite(j,bonds[z])
                    bond1::Int=bonds[z].site1.num;
                    bond2::Int=bonds[z].site2.num;
                    a=bonds[z].isNear ? J : J2;
                    #flip the bits at those relevant places
                    b::Int=flipBits(bond1-1, bond2-1,i);
                    push!(cols, b+1);
                    push!(rows, i+1);
                    push!(values, (1/2)*a/4);
                    if(eigmethod=="full")
                        H[i+1,b+1]=a/4;
                    end
                end
            end

        end
        #println("i count ", i, count);
    end
    if(eigmethod!="full")
        H=sparse(rows, cols, values);
    end
    if(N<=1)
        println("num sites, ", N, ", H", Matrix(H))
    end
    return H;
end

#TODO: update this and update the other one
#field is in the x direction!!! USE THIS ONE
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
        for j=1:N
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
        for j=1:N
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
        for j=1:N
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
        for j=1:N
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
        for j=1:N
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


#hamiltonian is sxisxj + h sz (field is in z direction)
function constructSusceptibilityHamiltonianFailed(states, bonds, N, J, J2, eigmethod, randList, randList2, sites2)
    #loop through ALL the possible states...
    println("bis");
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];

    list=states[1];
    map=states[2];
    #println(states);
    local H;
    if(eigmethod=="full")
        H=zeros(Float64, length(list),length(list));
    end
    #H::Matrix{Float64}=zeros(Int, length(list),length(list));

    for i=1:length(list)
        #NOW loop through all the possible SITES
        #flip bit
        #find state with that flipped bit

        count=0;
        for j=1:N
            #for this particular site, find ALL of its bonds
            #println("bonds length", length(bonds));
            theThing= getTi(j-1,list[i]) == 1 ? randList[j]/2 : -randList[j]/2;
            #println((1/2)*theThing);
            push!(rows, i);
            push!(cols, i);
            push!(values, theThing);
            if(eigmethod=="full")
                H[ilist[i],list[i]+1]+=theThing;
            end
            #perturbative field, sites2 is a list of sites to add the field in
                b::Int=flipBit(j-1, list[i]);
                t::Int=map[b];
                push!(rows, i);
                push!(cols, t);
                push!(values, ((-1)^(sites2[j]))*randList2[j]/2);
                if(eigmethod=="full")
                    H[i,t]=randList2[j]/2;
                end
            #perturbative field code end
            #bonds
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
                    b=flipBits(bond1-1, bond2-1,list[i]);
                    t=map[b];
                    a= bonds[z].isNear ? J : J2;
                    if(eigmethod=="full")
                        H[i,t]=a/4;
                    end
                    push!(rows, i);
                    push!(cols, t);
                    push!(values, (1/2)*a/4);
                end
            end

        end
        #println("i count ", i, count);
    end
    if(eigmethod!="full")
        H=sparse(rows, cols, values);
    end
    println("bisbhsibs");
    println(H);
    return H;
end


#states, bonds, N, J, J2, eigmethod, randList, randList2, sites2
function constructSusceptibilityHamiltonian(bonds, N, J, J2, eigmethod, randList, randList2, sites2)
    #loop through ALL the possible states...
    #println(states);
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];
    local H;

    if(eigmethod=="full")
        H=zeros(Float64, 2^(N), 2^(N));
    end
    #H::Matrix{Float64}=zeros(Int, length(list),length(list));

    for i=0:2^(N)-1
        #NOW loop through all the possible SITES
        count=0;
        for j=1:N
            #for this particular site, find ALL of its bonds
            #println("bonds length", length(bonds));
            theThing= getTi(j-1,i) == 1 ? randList[j]/2 : -randList[j]/2;
            #1/2 for double counting
            #println((1/2)*theThing);
            push!(cols, i+1);
            push!(rows, i+1);
            push!(values, theThing);
            if(eigmethod=="full")
                H[i+1,i+1]+= theThing;
            end
            #perturbative field
            b::Int=flipBit(j-1, i);
            push!(rows, i+1);
            push!(cols, b+1);
            push!(values, ((-1)^(sites2[j]))*randList2[j]/2);
            if(eigmethod=="full")
                H[i+1,t+1]=randList2[j]/2;
            end

            for z=1:length(bonds)
                #get the number of the thing it is bonding with!
                if containsSite(j,bonds[z])
                    bond1::Int=bonds[z].site1.num;
                    bond2::Int=bonds[z].site2.num;
                    a=bonds[z].isNear ? J : J2;
                    #flip the bits at those relevant places
                    b=flipBits(bond1-1, bond2-1,i);
                    push!(cols, b+1);
                    push!(rows, i+1);
                    push!(values, (1/2)*a/4);
                    if(eigmethod=="full")
                        H[i+1,b+1]=a/4;
                    end
                end
            end

        end
        #println("i count ", i, count);
    end
    if(eigmethod!="full")
        H=sparse(rows, cols, values);
    end
    if(N<=1)
        println("num sites, ", N, ", H", Matrix(H))
    end
    return H;
end

function constructTransverseHamiltonianNoSymmetrySxMeanField(bonds, N, J, J2, eigmethod, originalH, newFieldH)
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];

    if(eigmethod=="full")
        H=zeros(Float64, 2^(N), 2^(N));
    end
    #H::Matrix{Float64}=zeros(Int, length(list),length(list));
    for i=0:2^(N)-1
        #NOW loop through all the possible SITES
        for j=1:N
            #mean field theory field calculation
            theThing= getTi(j-1,i) == 1 ? newFieldH[j]/2 : -newFieldH[j]/2;
            #1/2 for double counting
            #println((1/2)*theThing);
            push!(cols, i+1);
            push!(rows, i+1);
            push!(values, theThing);
            if(eigmethod=="full")
                H[i+1,i+1]-= theThing;
            end
            #original field
            b::Int=flipBit(j-1, i);
            push!(rows, i+1);
            push!(cols, b+1);
            push!(values, -originalH[j]);
            if(eigmethod=="full")
                #originally just =
                H[i+1,b+1]=-originalH[j];
            end
            #for this particular site, find ALL of its bonds
            for z=1:length(bonds)
                #get the number of the thing it is bonding with!
                if containsSite(j,bonds[z])
                    bond1::Int=bonds[z].site1.num;
                    bond2::Int=bonds[z].site2.num;
                    local theThing;
                    a= bonds[z].isNear ? J : J2;
                    #println("bond1", bond1);
                    #println("bond2", bond2);

                    if(getTi(bond1-1,i)==getTi(bond2-1,i))
                        theThing=a/2;
                    else
                        theThing=-a/2;
                    end
                    push!(rows, i+1);
                    push!(cols, i+1);
                    push!(values, theThing);
                    if(eigmethod=="full")
                        H[i+1,i+1]+= theThing;
                    end
                end
            end
        end
    end
    if(eigmethod!="full")
        H=sparse(rows, cols, values);
    end
    if(N<=1)
        println("num sites, ", N, ", H", Matrix(H))
    end
    return H;
end


#input a d array of random on site values
function build_anderson_hamiltonian_1d(d, bonds, N, t)
    H=zeros(Float64, N, N)
    for i=1:N
        H[i, i]=d[i]
    end
    for i=1:length(bonds)
        bond1::Int=bonds[i].site1.num;
        bond2::Int=bonds[i].site2.num;
        H[bond1,bond2]=-t
        H[bond2,bond1]=-t
    end
    return H
end

# build the anderson model hamiltonian in number basis 
function build_anderson_hamiltonian_1d_interactions(d, bonds, N, t)
    H=zeros(Float64, N, N)
    for i=1:N
        H[i, i]=d[i]
    end
    for i=1:length(bonds)
        bond1::Int=bonds[i].site1.num;
        bond2::Int=bonds[i].site2.num;
        H[bond1,bond2]=-t
        H[bond2,bond1]=-t
    end
    return H
end

#hx and hz must be vectors

function constructTransverseHamiltonianNoSymmetrySxSz(bonds, N, J, J2, eigmethod, hx, hz)
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];

    if(eigmethod=="full")
        H=zeros(Float64, 2^(N), 2^(N));
    end
    #H::Matrix{Float64}=zeros(Int, length(list),length(list));
    for i=0:2^(N)-1
        #NOW loop through all the possible SITES
        for j=1:N
            b::Int=flipBit(j-1, i);
            push!(rows, i+1);
            push!(cols, b+1);
            push!(values, -hx[j]);

            push!(rows, i+1);
            push!(cols, i+1);
            spin=getTi(j-1, i)==1 ? 1 : -1
            push!(values, -hz[j]*spin)
            #println(values)

            if(eigmethod=="full")
                H[i+1,b+1]+=-hx[j];
                H[i+1,i+1]+=-hz[j]*spin;
            end
            #for this particular site, find ALL of its bonds
            #println("length:"*string(length(bonds)))
            for z=1:length(bonds)
                #get the number of the thing it is bonding with!
                if containsSite(j,bonds[z])
                    bond1::Int=bonds[z].site1.num;
                    bond2::Int=bonds[z].site2.num;
                    local theThing;
                    a= bonds[z].isNear ? J : J2;

                    if(getTi(bond1-1,i)==getTi(bond2-1,i))
                        theThing=a/2;
                    else
                        theThing=-a/2;
                    end
                    push!(rows, i+1);
                    push!(cols, i+1);
                    push!(values, theThing);
                    if(eigmethod=="full")
                        H[i+1,i+1]+= theThing;
                    end
                end
            end
        end
    end
    if(eigmethod!="full")
        H=sparse(rows, cols, values);
    end
    if(N<=1)
        println("num sites, ", N, ", H", Matrix(H))
    end
    return H;
end

# J1 (XX + YY) + J2 (ZZ) + hz Z (hz is a random field)
function constructTransverseHamiltonianNoSymmetrySxSzSy(bonds, N, J1, J2, eigmethod, hz)
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];

    if(eigmethod=="full")
        H=zeros(Float64, 2^(N), 2^(N));
    end
    #H::Matrix{Float64}=zeros(Int, length(list),length(list));
    for i=0:2^(N)-1
        #NOW loop through all the possible SITES
        for j=1:N
            # For the Z field
            push!(rows, i+1);
            push!(cols, i+1);
            spin=getTi(j-1, i)==1 ? 1 : -1
            push!(values, hz[j]*spin)
            #println(values)

            if(eigmethod=="full")
                H[i+1,i+1]+=hz[j]*spin;
            end
            #for this particular site, find ALL of its bonds
            for z=1:length(bonds)
                #get the number of the thing it is bonding with!
                if containsSite(j,bonds[z])
                    bond1::Int=bonds[z].site1.num;
                    bond2::Int=bonds[z].site2.num;
                    
                    # spin flipped at both sides (for the XX term)
                    flipped_state::Int=flipBit(bond2-1, flipBit(bond1-1, i));

                    HXX=J1/2
                    local HYY;
                    local HZZ;
                    # XX produces 1 no matter what
                    # up up: YY -1 , up down or down up: i * -i = 1, down down: 1

                    # divide by 2 for double counting!
                    if(getTi(bond1-1,i)==getTi(bond2-1,i))
                        HYY=-J1/2;
                        HZZ=J2/2
                    else
                        HYY=J1/2;
                        HZZ=-J2/2
                    end
                    push!(rows, i+1);
                    push!(cols, flipped_state+1);
                    push!(values, HXX+HYY);

                    push!(rows,i+1);
                    push!(cols, i+1);
                    push!(values, HZZ);

                    if(eigmethod=="full")
                        H[i+1,i+1]+= HZZ;
                        H[i+1, flipped_state+1] += (HXX+HYY)
                    end
                end
            end
        end
    end
    if(eigmethod!="full")
        H=sparse(rows, cols, values);
    end
    return H;
end


# Method to construct 1D Hamiltonian of spinless fermions with nearest neighbor interactions, hoppings, and onsite potentials h
# N is the number of sites
# J is the interaction energies (size N-1 if pbc=false, otherwise size N). J[j] is the interaction between j, j+1
# t: hoppings at each link, size N-1 if pbc=false, otherwise size N. t[j] is hopping between j and j+1
# computational basis
# pbc= true by default
function construct_disordered_interacting_hamiltonian_nearest_neighbor(N, J, h, t, pbc=true)
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];

    # Loop through all the states in the computational basis
    for i=0:2^(N)-1
        # Loop through all the sites
        # Keep track of the digaonal energy terms associated with this particular state 
        energy_of_i = 0
        for j=1:N
            # First we deal with interactions
            # get the occupation at the jth position for the state i
            next_site=0

            # occupation at current site
            occupation=getTi(j-1, i)

            # occupation at the neighboring site
            occupation_next=0
            # get its neighbors occupation. If not PBC and last site, there's no neighbor. If PBC, it's neighbors with the first site
            if(j<N)
                occupation_next=getTi(j, i)
            else
                # if periodic boundary conditions, we get the first site
                if(pbc)
                    occupation_next=getTi(0, i)               
                end
            end

            # the nearest neighbor interaction is the multiplication of the occupations and the interaction strength at link i J[i]
            energy_of_i += J[j]*occupation*occupation_next

            # We add onsite potentials
            energy_of_i += + h[j]

            # We now deal with hopping
            next_site = 0
            # flip the occupation at site j
            hopping_state = flipBit(j-1, i)
            # push the hopping energy to the off-diagonal element
            push!(rows, i+1)
            push!(cols, hopping_state+1)
            push!(values, t[j])

            # and the hermitian conjugate
            push!(rows, hopping_state+1)
            push!(cols, i+1)
            push!(values, t[j])

        end
        push!(rows, i+1);
        push!(cols, i+1);
        push!(values, energy_of_i)
    end

    H=sparse(rows, cols, values);

    return H;
end



# Method to construct 1D Hamiltonian of spinless hardcore bosons with nearest neighbor interactions, hoppings, and onsite potentials h
# N is the number of sites
# J is the interaction energies (size N-1 if pbc=false, otherwise size N). J[j] is the interaction between j, j+1
# t: hoppings at each link, size N-1 if pbc=false, otherwise size N. t[j] is hopping between j and j+1
# computational basis
# pbc= true by default
# need to return list of largest off diagonal elements
# build a sparse Hamiltonian
function construct_disordered_interacting_hamiltonian_nearest_neighbor(N, J, h, t, pbc=true, method="full")
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];

    # get the off-diagonal terms with the largest coefficient
    max_index = argmax(t)

    if(method=="full")
        H=zeros(Float64, 2^(N), 2^(N));
    end

    off_diags = zeros(Float64, N)

    # need to find all states that correspond to a particular hopping, and also store the largest of-diagonal element
    index_list = [[] for i=1:N]

    # Loop through all the states in the computational basis
    for i=0:2^(N)-1
        # Loop through all the sites
        # Keep track of the digaonal energy terms associated with this particular state 
        energy_of_i = 0
        for j=1:N
            # First we deal with interactions
            # get the occupation at the jth position for the state i
            next_site=0

            # occupation at current site
            occupation=getTi(j-1, i)

            # occupation at the neighboring site
            occupation_next=0
            # get its neighbors occupation. If not PBC and last site, there's no neighbor. If PBC, it's neighbors with the first site
            if(j<N)
                occupation_next=getTi(j, i)
            else
                # if periodic boundary conditions, we get the first site
                if(pbc)
                    occupation_next=getTi(0, i)               
                end
            end

            # the nearest neighbor interaction is the multiplication of the occupations and the interaction strength at link i J[i]
            energy_of_i += J[j]*occupation*occupation_next

            # We add onsite potentials
            energy_of_i += + h[j]*occupation

            # We now deal with hopping
            next_site = 0
            if(!pbc && j==N)
                continue
            end
            # flip the occupation at site j, then i
            if(occupation + occupation_next == 1)
                next_ind = j == N ? 0 : j
                hopping_state = flipBit(next_ind, flipBit(j-1, i))        
                # push the hopping energy to the off-diagonal element
                push!(rows, i+1)
                push!(cols, hopping_state+1)
                push!(values, -t[j]/2)

                push!(rows, hopping_state+1)
                push!(cols, i+1)
                push!(values, -t[j]/2)


                if(i<hopping_state)
                    push!(index_list[j], (i+1, hopping_state+1))
                end
                
                # and the hermitian conjugate
                if(method=="full")
                    H[hopping_state+1, i+1]-=t[j]/2;
                    H[i+1,hopping_state+1]-=t[j]/2;
                end
            end
        end
        push!(rows, i+1);
        push!(cols, i+1);
        push!(values, energy_of_i)
        if(method=="full")
            H[i+1, i+1]=energy_of_i;
        end
    end
    
    if(method=="sparse")
            H=sparse(rows, cols, values);
    end

    return (H, index_list);
end

# Method to construct 1D Hamiltonian of spinless hardcore bosons with nearest neighbor interactions, hoppings, and onsite potentials h
# N is the number of sites
# J is the interaction energies (size N-1 if pbc=false, otherwise size N). J[j] is the interaction between j, j+1
# t: hoppings at each link, size N-1 if pbc=false, otherwise size N. t[j] is hopping between j and j+1
# computational basis
# pbc= true by default
# need to return list of largest off diagonal elements
# build a sparse Hamiltonian
function construct_disordered_interacting_hamiltonian_nearest_neighbor_half_filling(N, J, h, t, map_half_filling_states, pbc=true, method="full")
    cols::Vector{Int}=Int[];
    rows::Vector{Int}=Int[];
    values::Vector{Float64}=Float64[];

    # get the off-diagonal terms with the largest coefficient
    max_index = argmax(t)
    all_states=collect(Base.keys(map_half_filling_states))
    #println(all_states)

    if(method=="full")
        H=zeros(Float64, (length(all_states)), (length(all_states)));
    end

    off_diags = zeros(Float64, N)

    # need to find all states that correspond to a particular hopping, and also store the largest of-diagonal element
    index_list = [[] for i=1:N]

    # Loop through all the computational basis states in the half-filling sector
    for i in all_states
        # Loop through all the sites
        # Keep track of the digaonal energy terms associated with this particular state 
        energy_of_i = 0
        i_index=map_half_filling_states[i]
        for j=1:N
            # First we deal with interactions
            # get the occupation at the jth position for the state i
            next_site=0

            # occupation at current site
            occupation=getTi(j-1, i)

            # occupation at the neighboring site
            occupation_next=0
            # get its neighbors occupation. If not PBC and last site, there's no neighbor. If PBC, it's neighbors with the first site
            if(j<N)
                occupation_next=getTi(j, i)
            else
                # if periodic boundary conditions, we get the first site
                if(pbc)
                    occupation_next=getTi(0, i)               
                end
            end

            # the nearest neighbor interaction is the multiplication of the occupations and the interaction strength at link i J[i]
            energy_of_i += J[j]*occupation*occupation_next

            # We add onsite potentials
            energy_of_i += + h[j]*occupation

            # We now deal with hopping
            next_site = 0

            if(!pbc && j==N)
                continue
            end
            # flip the occupation at site j, then i
            if(occupation + occupation_next == 1)
                # note that j is the next site
                next_ind = j == N ? 0 : j
                hopping_state = flipBit(next_ind, flipBit(j-1, i)) 

                hopping_state_index=map_half_filling_states[hopping_state]
                # push the hopping energy to the off-diagonal element
                push!(rows, i_index)
                push!(cols, hopping_state_index)
                push!(values, -t[j]/2)

                push!(rows, hopping_state_index)
                push!(cols, i_index)
                push!(values, -t[j]/2)

                if(i<hopping_state)
                    push!(index_list[j], (i_index, hopping_state_index))
                end
                
                # and the hermitian conjugate
                if(method=="full")
                    H[hopping_state_index, i_index]-=t[j]/2;
                    H[i_index,hopping_state_index]-=t[j]/2;
                end
            end
        end
        push!(rows, i_index);
        push!(cols, i_index);
        push!(values, energy_of_i)
        if(method=="full")
            H[i_index, i_index]=energy_of_i;
        end
    end
    
    if(method=="sparse")
            H=sparse(rows, cols, values);
    end

    return (H, index_list);
end