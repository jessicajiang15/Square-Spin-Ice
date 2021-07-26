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

function bondListTwo()
    list=bond[];
    site1=site(1, 0, 0);
    site2=site(2, 1, 0);
    site3=site(3, 0, 1);
    site4=site(4, 1, 1);

    bond1=bond(site1, site2);
    bond2=bond(site1, site3);
    bond3=bond(site2, site4);
    bond4=bond(site3, site4);
    bond5=bond(site1, site4);
    bond6=bond(site2, site3);
    push!(list, bond1);
    push!(list, bond2);
    push!(list, bond3);
    push!(list, bond4);
    push!(list, bond5);
    push!(list, bond6);

    return list;

end

function generateCheckerboardPlaquettes(N)
    plaquetteList=Int[];
    for i=1:N*N-N
        if(((i-1)%2!=0&&(i-1)÷N%2==0||(i-1)%2==0&&(i-1)÷N%2!=0)&&!isLast(i, N))
            push!(plaquetteList, i);
        end
    end
    return plaquetteList;
end

function isLast(i, N)
    return (i%N==0);
end

function generateCheckerboardNoCrossPlaquettes(N)
    plaquetteList=Int[];
    for i=1:(N*N)-N
        if((((i-1)%2==0&&(i-1)÷N%2==0)||((i-1)%2!=0&&(i-1)÷N%2!=0))&&!isLast(i, N))
            push!(plaquetteList, i);
        end
    end
    return plaquetteList;
end

function plaquetteIndicies(plaquetteIndex, N)
    return Int[plaquetteIndex, plaquetteIndex+1, plaquetteIndex+N, plaquetteIndex+N+1]
end

function getPlaquetteNumber(index, N)
    return getRow(index, N)%2==0 ? 1 : 0;
end


function getRow(index, N)
    return ((index-1)÷N)+1;
end


function generateListsofPlaquetteIndicies(N)
    plaquette=generateCheckerboardNoCrossPlaquettes(N);
    plaquetteLists=Vector{Int}[];
    for i=1:length(plaquette)
        push!(plaquetteLists, plaquetteIndicies(plaquette[i], N));
    end
    return plaquetteLists;
end
