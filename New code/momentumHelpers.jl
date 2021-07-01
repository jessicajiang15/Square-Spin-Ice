include("momentumBitOperations.jl")
include("necessaryBitOperations.jl")
include("necessaryStructures.jl")


function referenceStates(N)
    println("yes");
    mask=calculateMask(N);
    all::Array{Any}=Any[];
    #list of all reference states
    refStates::Array{refState}=refState[];

    #maps an integer to its corresponding state
    refMap::Dict{Int, state}=Dict{Int, state}();
    #for each reference state there are a set of compatible momenta
    #youll also need to know each state's periodicity
    ct=0;
    for i=0:2^(N*N)-1
                local temp=countPeriodicity(i, N, refMap);
                if(temp[3].periodicity==2)
                    ct+=1;
                end
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
        println("JANKNSJKANJSKA: ", ct);
        push!(all, refStates);
        push!(all, refMap);
        return all;
end

function containsSite(i, bond)
    return bond.site1.num==i||bond.site2.num==i;
end

function isViable(p, N, ref)
    return mod(p,(N/ref.periodicity))==0;
    #return mod(p, N*N÷ref.periodicity)==0;
    #return p%
end

function findViableRefStates(m, N, ref)
    viables::Array{refState}=refState[];
    all::Array{Any}=Any[];
    #links ref state to its numbering
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
    return 2*pi*psector/(N);
end


function getStateInfo(state, N, refMap)
    all::Array{Any}=Any[];
    isRefState::Bool=false;
    local ref::refState2d;
    yes::Bool=false;
    check::Int64=state;
    yShift::Int=0;
    xShift::Int=0;
    temp::Int=state;
    countUnique::Int=0;
    lol::Bool=false;
    for x=0:N-1
        xShift = x;
        for y=0:N-1
            if(x==0&&y==0)
                continue;
            end
            yShift = y;

            temp=rotateBits(x, y, state, N);
            #only count a state if youre a reference state– gives how many
            #unique states there are corresponding to particular ref state
            #once you get back to the ref state, and youre not one, all we need
            #is the translations that bring you to the ref state!
            if(temp==check&&yes)
                break;
            else
                countUnique+=1;
            end
                if(temp<state&&!yes)
                    ref=refMap[temp].ref;
                    isRefState=false;
                    check=ref.state;
                    yes=true;
                    lol=true;
                end
        end
    end
    if(!yes)
        ref=refState2d(state, xShift, yShift, countUnique);
        lol=true;
    end
    push!(all, isRefState);
    #shifts to get to ref state
    push!(all, xShift);
    push!(all, yShift);
    #corresponding ref state
    push!(all, ref);
    push!(all, countUnique);

end


function referenceStatesXY(N)
    println("yes");
    mask=calculateMask(N);
    all::Array{Any}=Any[];
    #list of all reference states
    refStates::Array{refState2d}=refState[];
    #maps an integer to its corresponding state
    refMap::Dict{Int, state2d}=Dict{Int, state2d}();
    #for each reference state there are a set of compatible momenta
    #youll also need to know each state's periodicity
    ct=0;
    for i=0:2^(N*N)-1
        local arr=getStateInfo(i, N, refMap);
        refMap[i]=state2d(i, arr[4], arr[2], arr[3]);
        if(arr[1])
            push!(refStates, arr[4]);
        end
    end
        push!(all, refStates);
        push!(all, refMap);
        return all;
end


function isViable2d(px, py, N, ref)
    theFactor::Complex{Int}=0;
    for i=0:N-1
        for j=0:N-1
            shifted=rotateBits(i, j, ref.state, N);
            if(shifted==ref)
                tf=exp(-1im*(i*calculateMomentum(px, N)+j*calculateMomentum(px, N)));
                if(tf<0)
                    return 0;
                end
                theFactor+=tf;
            end
        end
    end
        return theFactor;
    end



function getViableStates2d(px, py, N, ref)
    all::Array{Any}=Any[];
    viables::Array{refState2d}=refState2d[];
    dict::Dict{Int, Int}=Dict{Int, Int}();
    #state to its "sum"
    sumsOfPhaseFactors::Dict{Int, Complex{Int}}=Dict{Int, Complex{Int}}();
    count::Int=0;
    for i=1:length(ref)
        temp::Int=isViable2d(px, py, N, ref[i]);
        sumsOfPhaseFactors[i]=temp;
        if(temp!=0)
            count+=1;
            dict[ref[i].state]=count;
            push!(viables, ref[i]);
        end
    end
    push!(all, viables);
    push!(all, dict);
    push!(all, sumsOfPhaseFactors);
    return all;
end


function generateAllMomenta(N)
    momenta::Array{momentum}=momentum[];
    i=0;
    j=0;
    while(i<=N/2)
        while(j<=N/2)
            local p::momentum=momentum(i, j);
            push!(momenta, p);
            j+=1;
        end
        i+=1;
    end
    return momenta;
end
