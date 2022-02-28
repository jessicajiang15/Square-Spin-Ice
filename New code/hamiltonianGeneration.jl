include("bondGeneration.jl");
include("momentumHelpers.jl")
include("reflectionHelper.jl")
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
                H[i+1,i+1]+= theThing;
            end
            #original field
            b::Int=flipBit(j-1, i);
            push!(rows, i+1);
            push!(cols, b+1);
            push!(values, originalH[j]/2);
            if(eigmethod=="full")
                #originally just =
                H[i+1,b+1]=originalH[j]/2;
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
