include("reflectionBitOperations.jl")

function setState(refState, what, map)
    refState=map[what].refState;
    #println(refState);
    if(what==refState.state)
        return true;
    end
    return false;
end

function getReflectionStateInfo(state, N, map)
    temp::Array{Any}=Any[];
    isRefState::Bool=true;
    x::Int=0;
    y::Int=0;
    count::Int=1;

    refState=0;
    st=0;
    refX::Int=0;
    refY::Int=0;
    refXY::Int=0;
    xC::Int=0;
    yC::Int=0;
    xyC::Int=0;


    refX=reflectX(state, N);
    #println(refX);
    xC = (refX==state) ? 1 : -1;
    count= xC==-1 ? count+1 : count;
    if(refX<state)
        isRefState=false;
        refState=map[refX].refState;
        if(refState.state==refX)
            x=1;
            y=0;
            #println("WTF",refState);
            theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y,isRefState);
            return temp;
        end

    end
    #println(refX);

    refY=reflectY(state, N);
    yC =refY==state ? 1 : -1;
    count= yC==-1&&refY!=refX ? count+1 : count;

    if(refY<state&&isRefState)
            isRefState=false;
            refState=map[refY].refState;
            if(refY==refState.state)
                x=0;
                y=1;
                theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y,isRefState);
                return temp;
            end
        end
    if(!isRefState&&refY==refState.state)
                x=0;
                y=1;
                theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y,isRefState);
                return temp;
        end


        refX=reflectX(state, N);
    refXY =reflectY(refX, N);

    xyC =refXY==state ? 1 : -1;
    count= xyC==-1&&refXY!=refX&&refXY!=refY ? count+1 : count;
    if(refXY<state)
        isRefState=false;
        refState=map[refXY].refState;
        if(refState.state==refXY)
            x=1;
            y=0;
            theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y,isRefState);
            return temp;
            end
        end
    if(!isRefState&&refXY==refState.state)
            x=1;
            y=1;
            theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y,isRefState);
            return temp;
        end
        theEnd(temp, st, refState, state, xC, yC, xyC, count, x, y,isRefState);

        #println(temp);

        if(!isRefState)
        println("lmao wut");
    end
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
    tot=0;
    for i=0:2^(N*N)-1
        stateInfo=getReflectionStateInfo(i, N, map);
        if(stateInfo[1])
            if(stateInfo[3].xC*stateInfo[3].yC!=stateInfo[3].xyC)
                println("bish ", i);
                tot+=1;
            end
            push!(refs, stateInfo[3]);
        end
        map[i]=stateInfo[2];
    end
    println("tot", tot);
    push!(all, refs);
    push!(all, map);
    return all;
end


function isViableReflection(lambda, ref)
    factor=1;
    factor= ref.xC==1 ? factor+lambda.lx : factor;
    factor= ref.yC==1 ? factor+lambda.ly : factor;
    factor= ref.xyC==1 ? factor+lambda.lx*lambda.ly : factor;
    return factor!=0;
    #return (1+lambda.lx*ref.xC+lambda.ly*ref.yC+lambda.ly*lambda.lx*ref.xC*ref.yC!=0)
    #return (lambda.lx==ref.xC&&lambda.ly==ref.yC&&lambda.ly*lambda.lx==ref.xC*ref.yC);
end

function findViableRefStatesReflection(refs, lambda, N)
    all=Any[];
    map::Dict{Int, Int}=Dict{Int, Int}();
    viables::Array{reflectionRefState}=reflectionState[];
    count=0;
    for i=1:length(refs)
        if(isViableReflection(lambda, refs[i]))
            count+=1;
            push!(viables, refs[i]);
            map[refs[i].state]=count;
        end
    end
    push!(all, viables);
    push!(all, map);
end


function generateAllReflections()
    reflections=lambda[];
    lambda1=lambda(1, -1);
    #for this to be true reflection in x is same, not in y tho, so 2 unique states, xy reflection must be a different state
    lambda2=lambda(-1, -1);
    #none of them are the same... so 3 unique states? but ig x and y could result in same... xy reflection must go back
    lambda3=lambda(-1, 1);
    lambda4=lambda(1, 1);
    push!(reflections, lambda1);
    push!(reflections, lambda2);
    push!(reflections, lambda3);
    push!(reflections, lambda4);
    return reflections;


end
