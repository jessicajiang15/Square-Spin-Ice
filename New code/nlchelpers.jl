using DelimitedFiles;
include("calculations and plotting.jl")

function getCluster(order, N)
    list=Int[];
    return list;
end

function getAllSubClusters(order, N)
    subClusterCount=Int[];
    for i=1:order-1
        count=countCluster(order, i, N);
        push!(subClusterCount, count);
    end
    return subClusterCount;
end

function getWeight(order, N)
    if(order==1)
        return 1;
    end
    property=0;
    total=0;
    #subclusterCounts=getAllSubClusters(order, N);
    for i=1:order-1
            total+=getWeight(i);
    end
    return property-total;
end

function getCount(order, N)
end


function allGraphProperties(graphs)
    for i=1:length(graphs)
    end
end

function calculateNumSquares(graph)

end

function visitNodes(startingNode)
    #returns number of nodes visited from going through the whole chain
    list=Int[];
    count=0;
    startingNode.visited=true;
    push!(list, startingNode.num);
    count+=1;
    for i=1:length(startingNode.neighbors)
        if(!startingNode.neighbors[i].visited)
            push!(list, startingNode.neighbors[i]);
            append!(list, visitNodes(startingNode.neighbors[i]));
        end
    end
    return list;
end

function resetVisited(graph)
    for i=1:length(graph.nodes)
        graph.nodes[i].visited=false;
    end
end

function isSquare(graph)
    return length(visitNodes(graph.nodes[1]))==4
end

function numSquares(graph)
    count=0;
    visited=Int[];
    for i=1:length(graph.nodes)
        if(! graph.nodes[i].num in visited)
            temp=visitNodes(graph.nodes[i]);
            if(length(temp)==4)
                count+=1;
            end
            append!(visited, temp);
        end
    end
    return count;
end



function calculateNumSubgraphs(superGraph, subGraph)

end

function readFromGraphFile()
    r=readdlm("graphs.txt"; skipblanks=false);
    graphs::Vector{graph}=graph[];
    i=1;
    println("size", size(r)[1]);
    while(i<=size(r)[1])
        push!(graphs, readNumGraphs(r, i));
        i+=9;
    end
    return graphs;
end

function readNumGraphs(r, startingIndex)
    count=1;
    graphNum=0;
    numSites=0;
    numNearestNeighborBonds=0;
    numFarNeighborBonds=0;
    numSquares=0;
    numPlaquettes=0;
    numSubgraphs=0;
    latticeConstant=0;
    nearBonds::Vector{bond}=bond[];
    farBonds::Vector{bond}=bond[];
    squares::Vector{square}=square[];
    plaquettes::Vector{square}=square[];
    siteNumbering::Vector{Int}=Int[];
    subgraphList::Vector{Int}=Int[];
    row=1;
    for i=startingIndex:startingIndex+6
        if(i>size(r)[1])
            break;
        end
        count=1;
            if(row==1)
                num=r[i, count];

                while(num!=""&&count<=length(r[i,:]))
                    #println("wtf", num);
                    if(count==1)
                        if(i==10)
                            println("num", num);
                        end
                        graphNum=num;
                    elseif(count==2)
                        numSites=num;
                    elseif(count==3)
                        numNearestNeighborBonds=num;
                    elseif(count==4)
                        numFarNeighborBonds=num;
                    elseif(count==5)
                        numSquares=num;
                    elseif(count==6)
                        numPlaquettes=num;
                    elseif(count==7)
                        numSubgraphs=num;
                    elseif(count==8)
                        latticeConstant=num/2;
                    end
                    #println("count", num);
                    count+=1;
                    if(count<=length(r[i,:]))
                        num=r[i, count];
                    end
                end
            elseif(row==2)
                    num=r[i, count];
                    which=1;
                    num1=0;
                    num2=0;
                    while(num!=""&&count<=length(r[i,:]))
                        if(which==1)
                            num1=num;
                        else
                            num2=num;
                            push!(nearBonds, bond(site(num1, 0, 0), site(num2, 0, 0)));
                        end
                        which= which == 2 ? 1 : 2;
                        count+=1;
                        if(count<=length(r[i,:]))
                            num=r[i, count];
                        end
                    end
            elseif(row==3)
                num=r[i, count];
                which=1;
                num1=0;
                num2=0;
                while(num!=""&&count<=length(r[i,:]))
                    if(which==1)
                        num1=num;
                    else
                        num2=num;
                        push!(farBonds, bond(site(num1, 0, 0), site(num2, 0, 0)));
                    end
                    which= which == 2 ? 1 : 2;
                    count+=1;
                    if(count<=length(r[i,:]))
                        num=r[i, count];
                    end
                end
            elseif(row==4)
                num=r[i, count];
                which=1;
                num1=0;
                num2=0;
                num3=0;
                num4=0;
                while(num!=""&&count<=length(r[i,:]))
                    if(which==1)
                        num1=num;
                    elseif(which==2)
                        num2=num;
                    elseif(which==3)
                        num3=num;
                    elseif(which==4)
                        num4=num;
                        push!(squares, square(site(num1, 0, 0), site(num2, 0, 0), site(num3, 0, 0), site(num4, 0, 0)));
                        which=0;
                    end
                    which+=1;
                    count+=1;
                    if(count<=length(r[i,:]))
                        num=r[i, count];
                    end
                end
            elseif(row==5)
                num=r[i, count];
                #println("the row ",row, ", ", r[i, :]);
                which=1;
                num1=0;
                num2=0;
                num3=0;
                num4=0;
                while(num!=""&&count<=length(r[i,:]))
                    if(which==1)
                        num1=num;
                    elseif(which==2)
                        num2=num;
                    elseif(which==3)
                        num3=num;
                    elseif(which==4)
                        num4=num;
                        push!(plaquettes, square(site(num1, 0, 0), site(num2, 0, 0), site(num3, 0, 0), site(num4, 0, 0)));
                        which=0;
                    end
                    which+=1;
                    count+=1;
                    if(count<=length(r[i,:]))
                        num=r[i, count];
                    end
                end
            elseif(row==6)
                index=findfirst(isequal(""), r[i, :]);
                index = index==nothing ? length(r[i,:]) : index;
                siteNumbering=r[i, 1:index-1];
            elseif(row==7)
                index=findfirst(isequal(""), r[i, :]);
                index = index==nothing ? length(r[i,:]) : index;
                subgraphList=r[i, 1:index-1];
            end
            row+=1;
    end
    return graph(graphNum, numSites, numNearestNeighborBonds, numFarNeighborBonds, numSquares, numPlaquettes, numSubgraphs, latticeConstant, nearBonds, farBonds, squares, plaquettes, siteNumbering, subgraphList);
end

#=
function readNumGraphs(file, numGraphs)
    lines=readlines(file);
    num=0;
    numNearestNeighborBonds=0;
    numFarNeighborBonds=0;
    numSquares=0;
    numPlaquettes=0;
    numSubgraphs=0;
    latticeConstant=0;
    bonds::Vector{bond}=bond[];
    squares::Vector{square}=square[];
    plaquettes::Vector{square}=square[];
    siteNumbering::Vector{Int}=Int[];
    subgraphList=Any[];
    #true if you JUST read a " ", false otherwise
    yes=false;
    for j=1:length(lines)
        str=lines[j]
            if(i==1)
                    number="";
                    count=1;
                    for i=firstindex(str):lastindex(str)
                        if(str[i]==" ")
                            if(!yes)
                                num=parse(Int64, number);
                            else
                                num=0;
                            end
                            if(count==1)

                            elseif(count==2)
                            elseif(count==3)
                            elseif(count==4)
                            elseif(count==5)
                            elseif(count==6)
                            elseif(count==7)
                            elseif(count==8)

                            end
                            count+=1;
                        else
                            number*=str[i];
                        end
                    end
                end
            elseif(i==2)

            elseif(i==3)

            elseif(i==4)

            elseif(i==5)

            elseif(i==6)
            end
        end
    #initialize graph here
    graph=0;
    return graph;
end

=#

function calculateBaseWeightSz(J, h, width)
    bonds::Vector{bond}=bond[];
    temp=calculateEigensystemTransverse(1, J, h, bonds,"lanczos", "one", h, width);
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    sz=calculateSz(eigensystem[2], temp[3][eigensystem[3]], 1);
end

function calculateBasePropertySz()

end

function getAllWeightsSz(num, graphs, J, h, width)
    weights::Vector{Float64}=Float64[];
    push!(weights, calculateBaseWeightSz(J, h, width))
    for i=1:num
        println("order: ", i);
        push!(weights, calculateWeightSz(i, graphs, weights, J, h, width));
    end
    return weights;
end
#max order 56
function calculateInfiniteLatticeSz(order, J, h, graphs, width)
    @time begin
    println("weights starting!!");
    @time begin
        weights=getAllWeightsSz(order, graphs, J, h, width);
    end
    sum=0;

    println("weights done!!! ");
    for i=1:order
        println("order: ", order);
        sum+=getAllWeights[i]*graphs[i].latticeConstant;
    end
end
    return sum;
end

function calculateWeightSz(num, graphs::Vector{graph}, weights, J, h, width)
    sum=0;
    theGraph=graphs[num];
    list=theGraph.subgraphList;
    bonds=append!(theGraph.nearBonds, theGraph.farBonds);
    temp=calculateEigensystemTransverse(theGraph.numSites, J, h, bonds,"lanczos", "one", h, width);
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    sz=calculateSz(eigensystem[2], temp[3][eigensystem[3]], theGraph.numSites);
    for i=1:length(list)
        sum+=weights[list[i]];
    end
    return sz-sum;

end
