include("momentumBitOperations.jl")
include("necessaryBitOperations.jl")
include("necessaryStructures.jl")

function isHermitian(M)
    count=0;
    println("rows", size(M)[1]);
    for i=1:size(M)[1]
        for j=1+count:size(M)[2]
            thing=M[i, j]-conj(M[j, i]);
            #println(thing);
            if(abs(real(thing))>10^-10||abs(imag(thing))>10^-10)
                return false;
            end
        end
         count+=1;
    end
    return true;
end

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

function has(n, list)
    for i=1:length(list)
        if(list[i]==n)
            return true;
        end
    end
    return false;
end


function getStateInfo(state, N, refMap)
    all::Array{Any}=Any[];
    isRefState::Bool=true;
    local ref::refState2d;
    yes::Bool=false;
    check::Int64=state;
    yShift::Int=0;
    xShift::Int=0;
    temp::Int=state;
    countUnique::Int=1;
    list=Int64[];
    push!(list, state);

    lol::Bool=false;
    #println("state: ", state);

    for x=0:N-1
        if(yes)
            break;
        end
        xShift = x;
        for y=0:N-1
            if(x==0&&y==0)
                continue;
            end
            yShift = y;

            temp=rotateBits(x, y, state, N);
            #println("temp, ", temp);
            #only count a state if youre a reference state– gives how many
            #unique states there are corresponding to particular ref state
            #once you get back to the ref state, and youre not one, all we need
            #is the translations that bring you to the ref state!
            if(temp==check&&!isRefState)
                yes=true;
                break;
            end
            if(isRefState&&(!has(temp, list)))
                countUnique+=1;
                push!(list, temp);
            end
                if(temp<state&&isRefState)
                    ref=refMap[temp].ref;
                    isRefState=false;
                    if(ref.state==temp)
                        yes=true;
                        break;
                    end
                    check=ref.state;
                end
        end
    end
    if(isRefState)
        xShift=0;
        yShift=0;
        ref=refState2d(state, xShift, yShift, countUnique);
        lol=true;
    end

    push!(all, isRefState);
    #shifts to get to ref state
    push!(all, xShift);
    push!(all, yShift);
    #corresponding ref state
    push!(all, ref);
    #push!(all, countUnique);

end


function referenceStatesXY(N)
    println("yes");
    count::Int=0;
    all::Array{Any}=Any[];
    #list of all reference states
    refStates::Array{refState2d}=refState[];
    #maps an integer to its corresponding state
    refMap::Dict{Int, state2d}=Dict{Int, state2d}();
    #for each reference state there are a set of compatible momenta
    ct=0;
    for i=0:2^(N*N)-1
        local arr=getStateInfo(i, N, refMap);
        #red x y
        refMap[i]=state2d(i, arr[4], arr[2], arr[3]);
        if(arr[1])
            #println("unqieu statds", arr[4].numUniqueSt);
            push!(refStates, arr[4]);
        end
    end
    println("total !!!! ", count);
    println("how many ref states: ", length(refStates));
        push!(all, refStates);
        push!(all, refMap);
        return all;
end


function isViable2d(px, py, N, ref)
    theFactor::Complex{Float64}=0;
    for i=0:N-1
        for j=0:N-1
            shifted=rotateBits(i, j, ref.state, N);
            if(shifted==ref.state)
                tf=exp(-1im*(i*calculateMomentum(px, N)+j*calculateMomentum(py, N)));
                if(ref.state==6)
                    println(i," ", j);
                    println(tf);
                end
                theFactor+=tf;
            end
        end
    end
        return theFactor;
    end



function getViableStates2d(px, py, N, ref)
    all::Array{Any}=Any[];
    ct=0;
    viables::Array{refState2d}=refState2d[];
    dict::Dict{Int, Int}=Dict{Int, Int}();
    #state to its "sum"
    sumsOfPhaseFactors::Dict{Int, Complex{Float64}}=Dict{Int, Complex{Float64}}();
    count::Int=0;
    println("LENGTH OF REF", length(ref));
    for i=1:length(ref)
        temp::Complex{Float64}=isViable2d(px, py, N, ref[i]);
        #println("phase sum", abs(real(temp)));
        if(abs(real(temp))>10^-10)
            count+=1;
            dict[ref[i].state]=count;
            sumsOfPhaseFactors[ref[i].state]=real(temp);
            push!(viables, ref[i]);
        end
    end
    push!(all, viables);
    push!(all, dict);
    push!(all, sumsOfPhaseFactors);
    println("ct:",ct);
    return all;
end


function generateAllMomenta(N)
    momenta::Array{momentum}=momentum[];
    i=0;
    j=0;
    while(i<=N-1)
        j=0;
        while(j<=N-1)
            push!(momenta, momentum(i, j));
            j+=1;
        end
        i+=1;
    end
    println("length momentum", length(momenta));
    return momenta;
end
