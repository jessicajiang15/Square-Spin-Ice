#circularly shifts i by n digits, for a binary digit with N digits
function circularBitShift(n, i, N)
    mask::Int=1<<N-1;
    return ((i << mod(n,(N)))|(i >> (N-mod(n,(N)))))&mask;

end

#INCLUDING end digit
function circularBitShiftSector(n, i, N, startDigit, endDigit)
    mask::Int = ~(~0 << ((endDigit-1) - (startDigit-1) + 1));
    shifted::Int = (i >> (startDigit-1)) & mask;
    circShift::Int=circularBitShift(n, shifted, endDigit-startDigit+1);#rotate these digits
    circShift=circShift<<(startDigit-1);
    maskt=(((1 << (startDigit- 1)) - 1)⊻((1 << (endDigit)) - 1));
    return (i&(~(maskt)))|circShift;
end

function getRotatedBits(n, i, N, startDigit, endDigit)
    mask::Int = ~(~0 << ((endDigit-1) - (startDigit-1) + 1));
    shifted::Int = (i >> (startDigit-1)) & mask;
circShift::Int=circularBitShift(n, shifted, endDigit-startDigit+1);#rotate these digits

circShift=circShift<<(startDigit-1);
return circShift;
end

function rotateAllBitsX(n, k, N)
    test=Int[];
    for i=0:N-1
        push!(test,getRotatedBits(n, k, N, i*N+1, i*N+N));
    end
    result=0;
    for i=1:N
        result=result|test[i];
    end
    return result;
end

function rotateYBits(n, k, startD,N)
    theNum::Int=0;
    endD::Int=N*(N-1)+startD;
    mask::Int=0;
    i::Int=startD;
    final::Int=0;
    count::Int=0;
    while(i<=endD)
        #picks out nth digit and puts in right place
        #tmask is the ith digit
        local tMask::Int=((1<<(i-1))&k)>>(i-1-count);
        #then shift it back down and or
        theNum=tMask|theNum;
        i+=N;
        count+=1;
    end
    circ::Int=circularBitShift(n, theNum, N);

    j=0;
    z=startD;
    while(j<N)
        #picks out nth digit and puts in right place
        #tmask is the ith digit
        local tMask::Int=((1)&(circ>>(j)))<<(z-1);
        #then shift it back down and or
        final=tMask|final;
        j+=1;
        z+=N;
    end

    return final;
end



function rotateAllBitsY(n, k, N)
    test=Int[];
    for i=1:N
        #i is the column
        #start from i, end at i*N+i?
        push!(test, rotateYBits(n, k, i, N));
    end

    result=0;
    for i=1:N
        result=result|test[i];
    end
    return result;
end

function calculateMask(N)
    mask=0;
    temp=N*N-1;
    while(temp>=0)
        mask+=2^temp;
        temp-=1;
    end
    return mask;
end

function isSame(refState, i, N)
    for n=1:N
        if(circularBitShift(n, i, N)==refState)
            return n;
        end
    end
    return -1;
end

function countPeriodicity(b, N, map)
    temp::Array{Any}=Any[];
    yes::Bool=true;
    count::Int=1;
    shiftsNeeded::Int=0;
    local refSt::refState;
    checkPoint::Int=b;
    local refSt;
    alr=false;
    curr::Int=rotateAllBitsX(1, b, N);
    while(curr!=checkPoint)
        if(curr<b&&!alr)
            yes=false;
            refSt=map[curr].ref;
            if(checkPoint==refSt.state)
                break;
            end
            checkPoint=refSt.state;
            alr=true;
        end
        count+=1;
        curr=rotateAllBits(1, curr, N);
    end
    if(yes)
        refSt=refState(count, b);
    end
    shiftsNeeded=count;
    #is it a ref state
    push!(temp, yes);
    #shifts to get to ref state
    push!(temp, shiftsNeeded);
    #corresponding ref state
    push!(temp, refSt);
    #the periodicity
    push!(temp, count);
    return temp;
end



function countPeriodicity(b, N, map)
    temp::Array{Any}=Any[];
    yes::Bool=true;
    count::Int=1;
    shiftsNeeded::Int=0;
    local refSt::refState;
    checkPoint::Int=b;
    local refSt;
    curr::Int=rotateAllBitsX(1, b, N);
    while(curr!=checkPoint)
        if(curr<b)
            yes=false;
            refSt=map[curr].ref;
            if(checkPoint==refSt.state)
                break;
            end
            checkPoint=refSt.state;
        end
        count+=1;
        curr=rotateAllBitsX(1, curr, N);
    end
    if(yes)
        refSt=refState(count, b);
    end
    shiftsNeeded=count;
    #is it a ref state
    push!(temp, yes);
    #shifts to get to ref state
    push!(temp, shiftsNeeded);
    #corresponding ref state
    push!(temp, refSt);
    #the periodicity
    push!(temp, count);
    return temp;
end


function rotateBits(x, y, state, N)
    temp=rotateAllBitsX(x, state, N);
    return rotateAllBitsY(y, temp, N);
end
