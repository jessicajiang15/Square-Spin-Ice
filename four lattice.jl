#assume NxN lattice
#gets the ith binary digit (Starting from 0) of the decimal number I
using SparseArrays
using Arpack
using LinearAlgebra

struct site
    num::Int
    x::Int
    y::Int
end

struct bond
    site1::site
    site2::site
end

function getTi(i, I)
    return mod(I÷(2^(i)), 2);
end

#flips the bits at the digits i and j
#credit: https://stackoverflow.com/questions/18247126/how-to-flip-a-bit-at-a-specific-position-in-an-integer-in-any-language/18247246
function flipBits(i, j, n)
    temp::Int=n⊻(1<<i);
    temp::Int=temp⊻(1<<j);
    return temp;
end


#credit: semibran at https://github.com/semibran/wrap-around/blob/master/index.js
function wrap(n, m)
      return n >= 0 ? n % m : (n % m + m) % m;
end

#which: true: check i for periodic, false: check j for periodic
function addBond(i, j, which, N, bonds, site0)
    newi::Int=which=="i"||which=="ij" ? wrap(i, N) : i;
    newj::Int=which=="j"||which=="ij" ? wrap(j, N) : j;
    tempcount::Int=newj+1+(newi)*N;
    sit::site=site(tempcount, newi,newj);
    bnd::bond=bond(site0, sit);
    push!(bonds, bnd);
end


#assumes list is sorted and contains integers. finds the target a
function binarySearch(a,list, start, fin)
    if(start==fin)
        return nothing
    mid::Int=(fin÷2-start÷2);
    newstart::Int=start;
    newfin::Int=fin;
    cur::Intr=list[mid];
    if(curr==a)
        return mid
    else
        if(curr>a)
            newstart=mid;
        else
            newfin=mid;
        end
        return binarySearch(a, list, newstart, newfin);
    end

end

function containsSite(i, bond)
    return bond.site1.num==i||bond.site2.num==i;
end

function recalculateCount(count, i, j, inew, jnew, N)
    #i is row j is col
    tempcount=count;
    if(inew>i)
        tempcount+N;
    elseif(inew<i)
        tempcount-N;
    end
    if(tempcount<1)
        tempcount=N^2+tempcount;
    elseif(tempcount>N^2)
        tempcount=tempcount-N^2;
    end

    if(jnew>j)
        if(j%N==0)
            tempcount-=(N-1);
        else
            tempcount+=1;
        end
    elseif(jnew<j)
        if(j%N==1)
            tempcount+=(N-1);
        else
            tempcount-=1;
        end
    end
return tempcount;

end

function bondList(N)
    bonds::Array{bond}=bond[];
    if(N<=1)
        return bonds;
    end
    count::Int=0;
    for i=0:N-1
        for j=0:N-1
            count+=1;
            #bug: ened to recaluate count given current count and relative position
            site0::site=site(count, i, j);
            if(i!=N-1)
                addBond(i+1, j, "i", N, bonds, site0);
            end
            if(i==0)
                addBond(i-1, j, "i", N, bonds, site0);
            end

            if(j!=N-1)
                addBond(i, j+1, "j", N, bonds, site0);
            end

            if(j==0)
                addBond(i, j-1, "j", N, bonds, site0);
            end

            if(j!=N-1)
                addBond(i+1, j+1, "ij", N, bonds, site0);
            end

            if(j==0)
                addBond(i+1, j-1, "ij", N, bonds, site0);
            end

            if(j==0)
                addBond(i-1, j-1, "ij", N, bonds, site0);
            end
            if(j!=N-1)
                addBond(i-1, j+1, "ij", N, bonds, site0);
            end

        end
    end
    println(length(bonds));
    return bonds;
end

function generateHamiltonian()

end

#implementation of kernighan's algorithm (counting number of 1s)
function singleOutNUpSpins(N, range)
    println(N);
    temp::Array{Int}=Int[];
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
            push!(temp, i)
        end
    end
    return temp;
end

function calculateEigensystem(N)
    eigenvalues::Array{Any}=Any[];
    bonds::Array{bond}=bondList(N);
    for i=0:N*N
        if i==0
            continue;
        end
        spinUps::Array{Int}=singleOutNUpSpins(i, 2^(N*N));
        Htemp::SparseArrayCSC{Int}=constructHamiltonian(spinUps, bonds, N);
        #println(Htemp);
        eigtemp=eigs(Htemp, nev=length(spinUps));
        append!(eigenvalues, eigtemp[1]);
    end
    return eigenvalues;
end

function constructHamiltonian(states, bonds, N)
    #ok loop through ALL the possible states...
    println(length(states));
    H::SparseArrayCSC{Int}=spzeros(Int, length(states),length(states));
    for i=1:length(states)
        #NOW loop through all the possible SITES
        for j=1:N*N
            #for this particular site, find ALL of its bonds
            for j=1:length(bonds)
                #ok, we hit a bond, what do we do with it?
                #get the number of the thing it is bonding with!
                if containsSite(i,bonds[j])
                    bond1::Int=bonds[j].site1.num;
                    bond2::Int=bonds[j].site2.num;
                    #now, are the spins in those two places the same? if so, the
                    #diagonal entry at i,i is 1
                    if(getTi(bond1,i)==getTi(bond2, i))
                        H[i, i]+=1;
                    else
                        H[i, i]-=1;
                        #flip the bits at those relevant places

                        b::Int=flipBits(bond1, bond2,i);
                        t::Int=binarySearch(b, states, 1, length(states));
                        if t!=nothing
                            H[i,t]=2;
                        end
                    end

                end
            end
        end
            return H;
end


function runL()
    eigenvalues=calculateEigensystem(4);
    println(eigenvalues);
end

runL();
