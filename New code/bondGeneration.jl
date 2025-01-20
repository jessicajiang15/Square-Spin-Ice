include("necessaryBitOperations.jl");
include("necessaryStructures.jl");



#which: true: check i for periodic, false: check j for periodic
function addBond(i, j, which, N, bonds, site0, isNear)
    newi::Int=which=="i"||which=="ij" ? wrap(i, N) : i;
    newj::Int=which=="j"||which=="ij" ? wrap(j, N) : j;
    tempcount::Int=newj+1+(newi)*N;
    sit::site=site(tempcount, newi,newj);
    bnd::bond=bond(site0, sit, isNear);
    push!(bonds, bnd);
end

function addBondNoPBC(i, j, which, N, bonds, site0, isNear)
    if(i<0||j<0||i>(N-1)||j>(N-1))
        return;
    end
    tempcount::Int=j+1+(i)*N;
    sit::site=site(tempcount, i,j);
    bnd::bond=bond(site0, sit, isNear);
    push!(bonds, bnd);
end

function upperRight(i, j, N, bonds, site0)
    if(j!=N-1)
        addBond(i-1, j+1, "ij", N, bonds, site0, false);
    end
end

function bottomRight(i, j, N, bonds, site0)
    if(j!=N-1)
        addBond(i+1, j+1, "ij", N, bonds, site0, false);
    end
end

function bottomLeft(i, j, N, bonds, site0)
    if(j==0)
        addBond(i+1, j-1, "ij", N, bonds, site0, false);
    end
end

function upperLeft(i, j, N, bonds, site0)
    if(j==0)
        addBond(i-1, j-1, "ij", N, bonds, site0, false);
    end
end



function upperRightNoPBC(i, j, N, bonds, site0)
    if(j!=N-1)
        addBondNoPBC(i-1, j+1, "ij", N, bonds, site0, false);
    end
end

function bottomRightNoPBC(i, j, N, bonds, site0)
    if(j!=N-1)
        addBondNoPBC(i+1, j+1, "ij", N, bonds, site0, false);
    end
end

function bottomLeftNoPBC(i, j, N, bonds, site0)
    if(j==0)
        addBondNoPBC(i+1, j-1, "ij", N, bonds, site0, false);
    end
end

function upperLeftNoPBC(i, j, N, bonds, site0)
    if(j==0)
        addBondNoPBC(i-1, j-1, "ij", N, bonds, site0, false);
    end
end

function bondListFrustrated(N)
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
                addBond(i+1, j, "i", N, bonds, site0, true);
            end
            if(i==0)
                addBond(i-1, j, "i", N, bonds, site0, true);
            end

            if(j!=N-1)
                addBond(i, j+1, "j", N, bonds, site0, true);
            end

            if(j==0)
                addBond(i, j-1, "j", N, bonds, site0, true);
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
                addBond(i+1, j, "i", N, bonds, site0, true);
            end
            if(i==0)
                addBond(i-1, j, "i", N, bonds, site0, true);
            end

            if(j!=N-1)
                addBond(i, j+1, "j", N, bonds, site0, true);
            end

            if(j==0)
                addBond(i, j-1, "j", N, bonds, site0, true);
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
                addBond(i+1, j, "i", N, bonds, site0, true);
            end
            if(i==0)
                addBond(i-1, j, "i", N, bonds, site0, true);
            end

            if(j!=N-1)
                addBond(i, j+1, "j", N, bonds, site0, true);
            end

            if(j==0)
                addBond(i, j-1, "j", N, bonds, site0, true);
            end

            if(j!=N-1)
                addBond(i+1, j+1, "ij", N, bonds, site0, false);
            end

            if(j==0)
                addBond(i+1, j-1, "ij", N, bonds, site0, false);
            end

            if(j==0)
                addBond(i-1, j-1, "ij", N, bonds, site0, false);
            end
            if(j!=N-1)
                addBond(i-1, j+1, "ij", N, bonds, site0, false);
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

function generateCheckerboardNoCrossPlaquettesPlaq(N)
    plaquetteList=Int[];
    for i=1:(N*N-1)
        if((((i-1)%2!=0&&(i-1)÷N%2==0)||((i-1)%2==0&&(i-1)÷N%2!=0)))
            push!(plaquetteList, i);
        end
    end
    return plaquetteList;
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
    return Int[(getRow(plaquetteIndex, N)-1)*N+(plaquetteIndex-1)%4+1, (getRow(plaquetteIndex, N)-1)*N+(plaquetteIndex)%4+1, getRow(plaquetteIndex, N)==N ? (plaquetteIndex-1)%N+1 : plaquetteIndex+N, getRow(plaquetteIndex, N)==N ? (plaquetteIndex)%N+1 : plaquetteIndex+N+1]
end

function getPlaquetteNumber(index, N)
    if(getRow(index, N)%2==0)
        if(((index-1)%N)%2==0)
            return 0;
        else
            return 1;
        end
    else
        if(((index-1)%N)%2==0)
            return 1;
        else
            return 0;
        end
    end
end


function getPlaquetteNumberStart0(index, N)
    if(getRow(index, N)%2==0)
        if(((index-1)%N)%2==0)
            return 1;
        else
            return 0;
        end
    else
        if(((index-1)%N)%2==0)
            return 0;
        else
            return 1;
        end
    end
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

function generateListsofPlaquetteIndiciesFlip(N)
    plaquette=generateCheckerboardNoCrossPlaquettesPlaq(N);
    plaquetteLists=Vector{Int}[];
    for i=1:length(plaquette)
        push!(plaquetteLists, plaquetteIndicies(plaquette[i], N));
    end
    return plaquetteLists;
end

function getPlaquetteNumberList(N)
    nums=Int[];
    for i=1:N*N
        push!(nums, getPlaquetteNumber(i, N));
    end
    return nums;
end

function bonds1D(L, pbc)
    bonds::Array{bond}=bond[];
    for i=1:L-1
        sit0::site=site(i, i, 1)
        sit1::site=site(i+1, i+1, 1)
        bnd::bond=bond(sit0, sit1, true);
        push!(bonds, bnd);
    end
    if(pbc)
        sit0=site(L, L, 1)
        sit1=site(1, 1, 1)
        push!(bonds, bond(sit0, sit1,true)) 
    end
    return bonds
end



function bondListFrustratedNoPBC(N)
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
                    addBondNoPBC(i+1, j, "i", N, bonds, site0, true);
                end
                if(i==0)
                    addBondNoPBC(i-1, j, "i", N, bonds, site0, true);
                end

                if(j!=N-1)
                    addBondNoPBC(i, j+1, "j", N, bonds, site0, true);
                end

                if(j==0)
                    addBondNoPBC(i, j-1, "j", N, bonds, site0, true);
                end


                if(i%2==0)
                    #odd no \ and even no /
                    if(count%2==0)
                        upperRightNoPBC(i, j, N, bonds, site0);
                        bottomLeftNoPBC(i, j, N, bonds, site0);
                    else
                        upperLeftNoPBC(i, j, N, bonds, site0);
                        bottomRightNoPBC(i, j, N, bonds, site0);
                    end
                else
                    if(count%2!=0)
                        upperRightNoPBC(i, j, N, bonds, site0);
                        bottomLeftNoPBC(i, j, N, bonds, site0);
                    else
                        upperLeftNoPBC(i, j, N, bonds, site0);
                        bottomRightNoPBC(i, j, N, bonds, site0);
                    end
                end

        end
    end
    return bonds;
end

function isExternal(site, mapBondsNear, mapBondsFar)
    return 4-mapBondsNear[site]!=0||2-mapBondsFar[site]!=0;
end

#for now input N and try to get it to work!!
function calculateExternalFieldFactors(N, mapBondsNear, mapBondsFar)
    list=zeros(N*N);

    for i=1:N*N
        if(isExternal(i, mapBondsNear, mapBondsFar))
            num=getPlaquetteNumber(i, N);
            list[i]= num==0 ? -1 : 1;
        end
    end

    return list;
end
