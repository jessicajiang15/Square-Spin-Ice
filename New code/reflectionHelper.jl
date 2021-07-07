include("reflectionBitOperations.jl")

function setState(refState, what, map)
    refState=map[what].refState;
    println(refState);
    if(what==refState.state)
        return true;
    end
    return false;
end

function getReflectionStateInfo(state, N, map)
    temp::Array{Any}=Any[];
    local isRefState::Bool=true;
    local x::Int=0;
    local y::Int=0;
    local count::Int=1;

    local refState::reflectionRefState;
    local st::reflectionState;
    local refX::Int=0;
    local refY::Int=0;
    local refXY::Int=0;
    local xC::Int=0;
    local yC::Int=0;
    local xyC::Int=0;

    refX=reflectX(state, N);
    println(refX);
    xC = (refX==state) ? 1 : -1;
    #bug
    count= xC==1 ? count+1 : count;
    if(refX<state)
        isRefState=false;
        if(setState(refState, refX, map))
            x=1;
            y=0;
            println("WTF",refState);
            theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y,isRefState);
            return temp;
        end

    end
    println(refX);

    refY=reflectY(state, N);
    yC =refY==state ? 1 : -1;
    count= yC==1 ? count+1 : count;

    if(refY<state&&isRefState)
            isRefState=false;
            if(setState(refState, refX, map))
                x=0;
                y=1;
                theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y,isRefState);
                return temp;
            end
        elseif(!isRefState)
            if(refY==refState.state)
                x=0;
                y=1;
                theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y,isRefState);
                return temp;
            end
        end


    refXY =reflectY(refX, N);

    xyC =refXY==state ? 1 : -1;
    count= xyC==1 ? count+1 : count;
    if(refXY<state)
        isRefState=false;
        if(setState(refState, refX, map))
            x=1;
            y=0;
            theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y,isRefState);
            return temp;
            end
    elseif(!isRefState&&refY==refState.state)
            x=1;
            y=1;
            theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y,isRefState);
            return temp;
        end
        theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y,isRefState);

        println(temp);



    return temp;
end

function theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y, isRefState)
    if(isRefState)
        refState::reflectionRefState=reflectionRefState(state, xC, yC, xyC, count);
        st=reflectionState(state, refState, 0,0);
    else
        st=reflectionState(state, refState, x, y);
    end
    push!(temp, isRefState);
    push!(temp, st);
    push!(temp, refState);
end



function getReferenceStatesReflection(N)
    all::Array{Any}=Any[];
    map::Dict{Int, reflectionState}=Dict{Int, reflectionState}();
    refs::Array{reflectionRefState}=reflectionState[];
    for i=0:2^(N*N)-1
        stateInfo=getReflectionStateInfo(i, N, map);
        if(stateInfo[1])
            push!(refs, stateInfo[3]);
        end
        map[i]=stateInfo[2];
    end
    push!(all, refs);
    push!(all, map);
    return all;
end


function isViableReflection(lambda, ref)
    return (lambda.lx==ref.xC&&lambda.ly==ref.xY&&lambda.ly*lambda.lx==ref.xyC);
end

function findViableRefStatesReflection(refStates, lambda, N)
    all=Any[];
    map::Dict{Int, Int}=Dict{Int, Int}();
    viables::Array{reflectionRefState}=reflectionState[];
    for i=1:length(refStates)
        if(isViableReflection(lambda, refStates[i]))
            push!(viables, ref[i]);
            map[ref[i].state]=count;
        end
    end
    push!(all, viables);
    push!(all, map);
end


function generateAllReflections()
    reflections=lambda[];
    lambda1=lambda(1, -1);
    lambda2=lambda(-1, -1);
    lambda3=lambda(-1, 1);
    lambda4=lambda(1, -1);
    push!(reflections, lambda1);
    push!(reflections, lambda2);
    push!(reflections, lambda3);
    push!(reflections, lambda4);
    return reflections;


end
