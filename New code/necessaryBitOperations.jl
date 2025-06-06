using Distributions
using Random
function binarySearchIt(a, list)
    i=0;
    j=length(list)-1;
    mid=i+(j÷2-i÷2);
    while(j>=i)
        if(list[mid+1]==a)
            return mid+1;
        else
            if(a>list[mid+1])
                i=mid+1;
                mid=i+(j÷2-i÷2);
            else
                j=mid-1;
                mid=i+(j÷2-i÷2);
            end
        end
    end
    return -1;
end





function singleOutEvenOddSpins(isEven,range)
    println(range);
    #temp2 stores the map and the list
    temp2::Array{Any}=Any[];
    #temp stores the states
    temp::Array{Int}=Int[];
    #stores the states and their corresponding number
    temp1::Dict{Int,Int}=Dict{Int, Int}();
    in::Int=0;
    for i=0:range-1
        count::Int=0;
        n::Int=i;
        while n>0
            count+=1;
            n=n & (n-1);
        end
        if(isEven&&count%2==0)
            in+=1;
            temp1[i]=in;
            push!(temp, i)
        elseif(!isEven&&count%2!=0)
            in+=1;
            temp1[i]=in;
            push!(temp, i)
    end
end
push!(temp2, temp);
push!(temp2, temp1);
    return temp2;

end

function getTi(i, I)
    #println(2^i);
    return mod(I÷(2^(i)), 2);
end


#flips the bits at the digits i and j
#credit: https://stackoverflow.com/questions/18247126/how-to-flip-a-bit-at-a-specific-position-in-an-integer-in-any-language/18247246
function flipBits(i, j, n)
    temp::Int=n⊻(1<<i);
    temp=temp⊻(1<<j);
    return temp;
end

function flipBit(i, n)
    temp::Int=n⊻(1<<i);
    return temp;
end

#credit: semibran at https://github.com/semibran/wrap-around/blob/master/index.js
function wrap(n, m)
      return n >= 0 ? n % m : (n % m + m) % m;
end

function countBits(n)
    count::Int=0;
    temp::Int=n;
    while(temp>0)
        temp&=temp-1;
        count+=1;
    end
    return count;
end


function singleOutNUpSpins(N, range)
    #temp2 stores the map and the list
    temp2::Array{Any}=Any[];
    #temp stores the statess
    temp::Array{Int}=Int[];
    #stores the states and their corresponding number
    temp1::Dict{Int,Int}=Dict{Int, Int}();
    in::Int=0;
    for i=0:range-1
        count::Int=0;
        n::Int=i;
        while n>0
            count+=1;
            if count>N
                break;
            end
            n=n & (n-1);
        end
        if(count==N)
            in+=1;
            temp1[i]=in;
            push!(temp, i)

    end
end
println("THE STATES",length(temp));
push!(temp2, temp);
push!(temp2, temp1);
    return temp2;

end

function generateRandomh(hbar, width, bonds)
    temp=Normal(hbar, width);
    return length(bonds)==0 ? rand(temp, 1) : rand(temp, length(bonds));
end

#assume index 1 input
function getBit(i, n)
    return (n >> (i-1)) & 1
end