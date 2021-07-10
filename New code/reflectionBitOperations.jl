function reflectX(n, N)
    i=0;
    reset=0;
    currRow=1;
    curr=n;
    while(currRow<=N)
        for i=0:N÷2-1
             curr=swapBits(i+(currRow-1)*N, currRow*N-i-1, curr);
        end
        currRow+=1;
    end
    return curr;
end

function reflectY(n, N)
    i=0;
    reset=0;
    currCol=0;
    curr=n;
    while(currCol<N/2)
        for i=0:N-1
             curr=swapBits(currCol*N+i, (N*(N-1)+i)-currCol*N, curr);
        end
        currCol+=1;
    end
    return curr;
end


function swapBits(i, j, n)
    bit1=(n>>i)&1
    bit2=(n>>j)&1;
    x=bit1⊻bit2;
    x=(x<<i)|(x<<j);
    return n⊻x;
end
