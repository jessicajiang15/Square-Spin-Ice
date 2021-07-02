function reflectX(n, N)

end


function swapBits(i, j, n)
    bit1=(n>>i)&1
    bit2=(n>>j)&1;
    x=bit1⊻bit2;
    x=(x<<i)|(x<<j);
    return n⊻x;
end
