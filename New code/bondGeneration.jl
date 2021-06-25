include("necessaryBitOperations.jl");
include("necessaryStructures.jl");



#which: true: check i for periodic, false: check j for periodic
function addBond(i, j, which, N, bonds, site0)
    newi::Int=which=="i"||which=="ij" ? wrap(i, N) : i;
    newj::Int=which=="j"||which=="ij" ? wrap(j, N) : j;
    tempcount::Int=newj+1+(newi)*N;
    sit::site=site(tempcount, newi,newj);
    bnd::bond=bond(site0, sit);
    push!(bonds, bnd);
end

function upperRight(i, j, N, bonds, site0)
    if(j!=N-1)
        addBond(i-1, j+1, "ij", N, bonds, site0);
    end
end

function bottomRight(i, j, N, bonds, site0)
    if(j!=N-1)
        addBond(i+1, j+1, "ij", N, bonds, site0);
    end
end

function bottomLeft(i, j, N, bonds, site0)
    if(j==0)
        addBond(i+1, j-1, "ij", N, bonds, site0);
    end
end

function upperLeft(i, j, N, bonds, site0)
    if(j==0)
        addBond(i-1, j-1, "ij", N, bonds, site0);
    end
end

function bondListFrustrated(N)
    println("timing bondlist");
    @time begin
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

            if(i%2==0)
                #odd no \ and even no /
                if(count%2==0)
                    upperRight(i, j, N, bonds, site0);
                    bottomLeft(i, j, N, bonds, site0);
                else
                    upperLeft(i, j, N, bonds, site0);
                    bottomRight(i, j, N, bonds, site0);
                end
            else
                if(count%2!=0)
                    upperRight(i, j, N, bonds, site0);
                    bottomLeft(i, j, N, bonds, site0);
                else
                    upperLeft(i, j, N, bonds, site0);
                    bottomRight(i, j, N, bonds, site0);
                end
            end

        end
    end
    println(length(bonds));
end
    return bonds;
end


function bondListFourNeighbors(N)
    println("timing bondlist");
    @time begin
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
        end
    end
    println(length(bonds));
end
    return bonds;
end

function bondListEightNeighbors(N)
    println("timing bondlist");
    @time begin
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
end
    return bonds;
end
