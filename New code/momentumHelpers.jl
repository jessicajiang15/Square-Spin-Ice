include("momentumBitOperations.jl")
include("necessaryBitOperations.jl")
include("necessaryStructures.jl")

function getViableP(N, p)
    temp::Array{Float64}=Float64[];
    for m=0:p-1
        push!(temp, m/p);
    end
    return temp;
end

function checkCompatibility(ref, p, N)
    return mod(p, N*N÷ref.periodicity)!=0;
end

function referenceStates(N)
    mask=calculateMask(N);
    all::Array{Any}=Any[];
    #list of all reference states
    refStates::Array{refState}=refState[];

    #maps an integer to its corresponding state
    refMap::Dict{Int, state}=Dict{Int, state}();
    #for each reference state there are a set of compatible momenta
    #youll also need to know each state's periodicity

    for i=0:2^(N*N)-1
                local temp=countPeriodicity(i, N, refMap);
                if(temp[1])
                    #is a reference state
                    local tempSt::refState=refState(temp[4], i);
                    refMap[i]=state(i, tempSt, 0);
                    push!(refStates, tempSt);

                else
                    local shifts=temp[2];
                    refMap[i]=state(i, refMap[temp[3].state].ref, shifts);

                end
                #=
                local shifts=isSame(refStates[j].state,i, N);
                if(shifts!=-1)
                    refMap[i]=state(i, refStates[j], shifts);
                    continue;
                else
                    local temp::Int=countPeriodicity(i, N);
                    local arr::Array{Float64}=getViableP(i, N, temp);
                    local tempSt::refState=refState(temp, i, arr);
                    refMap[i]=state(i, tempSt, 0);
                    push!(refStates, tempSt);
                end

                if(length(refStates)<1)
                    local tempSt::refState=refState(1, i);
                    refMap[i]=state(i, tempSt, 0);
                    push!(refStates, tempSt);
                    continue;
                end
                =#
        end
        push!(all, refStates);
        push!(all, refMap);
        return all;
end

function containsSite(i, bond)
    return bond.site1.num==i||bond.site2.num==i;
end

function isViable(p, N, ref)
    return mod(p, N*N÷ref.periodicity)==0;
end

function findViableRefStates(m, N, ref)
    viables::Array{refState}=refState[];
    all::Array{Any}=Any[];
    dict::Dict{Int,Int}=Dict{Int, Int}();
    count::Int=0;
    for i=1:length(ref)
        if(isViable(m, N, ref[i]))
            count+=1;
            push!(viables, ref[i]);
            dict[ref[i].state]=count;
        end
    end
    push!(all, viables);
    push!(all, dict);
    return all;
end

function calculateMomentum(psector, N)
    return 2*pi*psector/(N*N);
end
