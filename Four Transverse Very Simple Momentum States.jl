using Arpack
using LinearAlgebra

#circularly shifts i by n digits, for a binary digit with N digits
function circularBitShift(n, i, N)
    return (i << n)|(i >> (N - n));
end

function isSame(refState, i, N)
    for n=0:N
        if(circularBitShift(n,i, N))
            return true;
        end
    end
    return false;
end

function referenceStates(N)
    referenceStates=Float64[];
    for i=0:2^N-1
            for j=0:length(referenceStates)
                if(isSame(referenceStates[j],i, N))
                    continue;
                end
            end
            push!(referenceStates, i);
        end
end
