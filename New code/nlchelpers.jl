using DelimitedFiles;
include("calculations and plotting.jl")

function getLastGraphNumOrder(order, graphs)
    temp=graphs[1].numSquares;
    count=1;
    while(count<length(graphs)&&temp<=order)
        count+=1
        temp=graphs[count].numSquares;
    end
    return count;
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
                            push!(nearBonds, bond(site(num1, 0, 0), site(num2, 0, 0), true));
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
                        push!(farBonds, bond(site(num1, 0, 0), site(num2, 0, 0), false));
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

function calculateBaseWeightSz(J, J2, h, width)
    bonds::Vector{bond}=bond[];
    temp=calculateEigensystemTransverse(1, J, J2, h, bonds,"lanczos", "one", h, width);
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    sz=calculateSz(eigensystem[2], temp[3][eigensystem[3]], 1);
    return sz;
end

function getAllWeightsSz(num, graphs, J, J2, h, width)
    local weights::Vector{Float64}=Float64[];
    base=calculateBaseWeightSz(J, J2, h, width);
    push!(weights, base)
    for i=1:num
        println("order: ", i);
        push!(weights, calculateWeightSz(i, graphs, weights, J, J2, h, width));
    end
    return weights;
end
#max order 56
function calculateInfiniteLatticeSz(order::Int, J, J2, h, graphs, width)
    @time begin
    println("weights starting!!");
    @time begin
        weights=getAllWeightsSz(order, graphs, J, J2, h, width);
    end
    sum=weights[1];

    println("weights done!!! ");
    for i=1:order
        #println("order: ", order);
        sum+=weights[i+1]*graphs[i].latticeConstant;
    end
end
    return sum;
end

function calculateWeightSz(num, graphs::Vector{graph}, weights, J, J2, h, width)
    sum=0;
    #the graph to calculate weight of
    theGraph=graphs[num];
    list=theGraph.subgraphList;
    temp=copy(theGraph.nearBonds);
    bonds=append!(temp, theGraph.farBonds);
    temp=calculateEigensystemTransverse(theGraph.numSites, J, J2, h, bonds,"lanczos", "one", h, width);
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    println("h: ", h);
    println("eigenvalues: ", eigenvalues);
    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    sz=theGraph.numSites*calculateSz(eigensystem[2], temp[3][eigensystem[3]], theGraph.numSites);
    for i=1:length(list)
        sum+=weights[list[i]+1];
    end
    sum+=theGraph.numSites*weights[1];
    return sz-sum;

end




function calculateWeightNoSubSz(num, graphs::Vector{graph}, J, J2, h, width)
    sum=0;
    #the graph to calculate weight of
    theGraph=graphs[num];
    list=theGraph.subgraphList;
    bonds=append!(theGraph.nearBonds, theGraph.farBonds);
    temp=calculateEigensystemTransverse(theGraph.numSites, J, J2, h, bonds,"lanczos", "one", h, width);
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    sz=theGraph.numSites*calculateSz(eigensystem[2], temp[3][eigensystem[3]], theGraph.numSites);
    return sz;
end

function getAllWeightsNoSubSz(num, graphs, J, J2, h, width)
    weights::Vector{Float64}=Float64[];
    base=calculateBaseWeightSz(J, J2, h, width);
    push!(weights, base)
    for i=1:num
        #println("order: ", i);
        push!(weights, calculateWeightNoSubSz(i, graphs, J, J2, h, width));
    end
    return weights;
end



function calculateBaseWeightEntanglement(J, J2, h, width)
    return 0;
end

function getAllWeightsEntanglement(num, graphs, J, J2, h, width)
    local weights::Vector{Float64}=Float64[];
    push!(weights, 0)
    #1st order
    push!(weights, 0);
    for i=2:num
        #println("order: ", i);
        push!(weights, calculateWeightEntanglement(i, graphs, weights, J, J2, h, width));
    end
    return weights;
end
#max order 56
function calculateInfiniteLatticeEntanglement(order::Int, J, J2, h, graphs, width)
    @time begin
    println("weights starting!!");
    @time begin
        weights=getAllWeightsEntanglement(order, graphs, J, J2, h, width);
    end

    println("weights done!!! ");
        sum=weights[1];
        for i=1:order
            #println("order: ", order);
            sum+=weights[i+1]*graphs[i].latticeConstant;
    end

end
    return sum;
end

function convertToList(squares)
    list=Vector{Int}[];
    for i=1:length(squares)
        temp=Int[];
        push!(temp, squares[i].num1.num);
        push!(temp, squares[i].num2.num);
        push!(temp, squares[i].num3.num);
        push!(temp, squares[i].num4.num);
        push!(list, sort(temp));
    end
    return list;
end

function calculateWeightEntanglement(num, graphs::Vector{graph}, weights, J, J2, h, width)
    sum=0;
    #the graph to calculate weight of
    theGraph=graphs[num];
    list=theGraph.subgraphList;

    squareList=convertToList(theGraph.squares);

    temp=copy(theGraph.nearBonds);

    bonds=append!(temp, theGraph.farBonds);
    temp=calculateEigensystemTransverse(theGraph.numSites, J, J2, h, bonds,"lanczos", "one", h, width);
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    entanglement=0;

    for i=1:length(squareList)
        listA=squareList[i];
        entanglement+=getEntanglementEntropy(eigensystem[2], temp[3][eigensystem[3]], listA, theGraph.numSites);
    end

    for i=1:length(list)
        sum+=weights[list[i]+1];
    end
    return entanglement-sum;
end




function calculateWeightNoSubEntanglement(num, graphs::Vector{graph}, J, J2, h, width)
    sum=0;
    #the graph to calculate weight of
    theGraph=graphs[num];
    list=theGraph.subgraphList;

    squareList=convertToList(theGraph.squares);

    temp=copy(theGraph.nearBonds);

    bonds=append!(temp, theGraph.farBonds);
    temp=calculateEigensystemTransverse(theGraph.numSites, J, J2, h, bonds,"lanczos", "one", h, width);
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    entanglement=0;

    for i=1:length(squareList)
        listA=(squareList[i]);
        entanglement+=getEntanglementEntropy(eigensystem[2], temp[3][eigensystem[3]], listA, theGraph.numSites);
    end
    return entanglement;

end

function getAllWeightsNoSubEntanglement(num, graphs, J, J2, h, width)
    weights::Vector{Float64}=Float64[];
    push!(weights, 0)
    push!(weights, 0);
    for i=2:num
        #println("order: ", i);
        push!(weights, calculateWeightNoSubEntanglement(i, graphs, J, J2, h, width));
    end
    return weights;
end



#HASJKJASJKANKSA


function calculateBaseWeightSpi(J, J2, h, width)
    bonds::Vector{bond}=bond[];
    temp=calculateEigensystemTransverseNoSymmetry(1, J, J2, h, bonds,"lanczos", "one", h, width, "H1");
    eigenvalue = temp[1][1]
    eigenvector = temp[2]
    sz=calculateSPiSzNew(eigenvector, temp[3][1], 1, Int[1]);
    return sz;
end

function getAllWeightsSpi(num, graphs, J, J2, h, width)
    local weights::Vector{Float64}=Float64[];
    base=calculateBaseWeightSpi(J, J2, h, width);
    push!(weights, base)
    for i=1:num
        #println("order: ", i);
        push!(weights, calculateWeightSpi(i, graphs, weights, J, J2, h, width));
    end
    return weights;
end
#max order 56

function calculateInfiniteLatticeSpi(orders::Vector{Int}, J, J2, h, graphs, width)
    @time begin
    println("weights starting!!");
    list=Float64[];
    @time begin
        weights=getAllWeightsSpi(orders[length(orders)], graphs, J, J2, h, width);
    end

    println("weights done!!! ");
    for order in orders
        sum=weights[1];
        for i=1:order
            #println("order: ", order);
            sum+=weights[i+1]*graphs[i].latticeConstant;
        end
        push!(list, sum);
    end

end
    return list;
end


function calculateWeightSpi(num, graphs::Vector{graph}, weights, J, J2, h, width)
    sum=0;
    #the graph to calculate weight of
    theGraph=graphs[num];
    list=theGraph.subgraphList;
    temp=copy(theGraph.nearBonds);
    bonds=append!(temp, theGraph.farBonds);
    temp=calculateEigensystemTransverseNoSymmetry(theGraph.numSites, J, J2, h, bonds,"lanczos", "one", h, width, "H1");
    eigenvector = temp[2]
    #eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    sz=theGraph.numSites*calculateSPiSzNew(eigenvector, temp[3][1], theGraph.numSites, theGraph.indicies);
    for i=1:length(list)
        sum+=weights[list[i]+1];
    end
    sum+=theGraph.numSites*weights[1];
    return sz-sum;
end




function calculateWeightNoSubSpi(num, graphs::Vector{graph}, J, J2, h, width)
    sum=0;
    #the graph to calculate weight of
    theGraph=graphs[num];
    list=theGraph.subgraphList;
    bonds=append!(theGraph.nearBonds, theGraph.farBonds);
    temp=calculateEigensystemTransverseNoSymmetry(theGraph.numSites, J, J2, h, bonds,"lanczos", "one", h, width, "H1");
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    sz=theGraph.numSites*calculateSPiSzNew(eigensystem[2], temp[3][eigensystem[3]], theGraph.numSites, theGraph.indicies);
    return sz;

end

function getAllWeightsNoSubSpi(num, graphs, J, J2, h, width)
    weights::Vector{Float64}=Float64[];
    base=calculateBaseWeightSpi(J, J2, h, width);
    push!(weights, base)
    for i=1:num
        #println("order: ", i);
        push!(weights, calculateWeightNoSubSpi(i, graphs, J, J2, h, width));
    end
    return weights;
end




function calculateInfiniteLatticeSz(orders::Vector{Int}, J, J2, h, graphs, width)
    @time begin
    println("weights starting!!");
    szs=Any[];
    @time begin
        weights=getAllWeightsSz(orders[length(orders)], graphs, J, J2, h, width);
        println("weights", weights);
    end

    println("weights done!!! ");
    for order in orders
        sum=weights[1];
        for i=1:order
            #println("order: ", order);
            sum+=weights[i+1]*graphs[i].latticeConstant;
        end
        push!(szs, sum);
    end

end
    return szs;
end

#here order=number of graphs
function calculateInfiniteLatticeSpi(order::Int, J, J2, h, graphs, width)
    @time begin
    println("weights starting!!");
    @time begin
        weights=getAllWeightsSpi(order, graphs, J, J2, h, width);
    end
    sum=weights[1];

    println("weights done!!! ");
    for i=1:order
        #println("order: ", order);
        sum+=weights[i+1]*graphs[i].latticeConstant;
    end
end
    return sum;
end






#lmao

function getBaseStatesList(order, graphs, hs, J, J2, width, hstep)
    hstep=(hs[length(hs)]-hs[1])/length(hs);
    ultList=Any[];
    listA=[];
        for h in hs
            theGraph=graphs[order];
            list=theGraph.subgraphList;
            temp=copy(theGraph.nearBonds);
            bonds=append!(temp, theGraph.farBonds);
            temp=calculateEigensystemTransverseNoSymmetry(1, J, J2, h, bonds,"lanczos", "one", h, width, "H1");
            eigenvalues = temp[1]
            eigenvectors = temp[2]
            local currState=Any[];
            push!(currState,eigenvectors)
            push!(currState, 0:1);
            push!(listA,currState)
        end
    push!(ultList, listA);
    return ultList;
end

#=
organizing fidelity:
1. need a list of all the "states" for any particular order (all the different hs)
2. loop through ALL those states and get a list of fidelities from that list
3. now, go onto the next order and do the same thing
=#

#returns a list of ALL the gs eigenvectors from hmin to hmax of the graph with num=order, with num many h values
#returns all states of differnet h at a particular order
function generateStatesList(order, graphs, hs, J, J2, width, hstep)
    states=[];
    vectors=[];
    ultimateList=[];
    for i=1:order
        listA=Any[];
        for h in hs
            theGraph=graphs[i];
            list=theGraph.subgraphList;
            temp=copy(theGraph.nearBonds);
            bonds=append!(temp, theGraph.farBonds);
            temp=calculateEigensystemTransverseNoSymmetry(theGraph.numSites, J, J2, h, bonds,"lanczos", "one", h, width, "H1");
            currState=Any[];
            push!(currState,temp[2])
            push!(currState, temp[3][1]);
            push!(listA,currState)
        end
        push!(ultimateList, listA);
    end
    return ultimateList;
end

#returns a list of ALL the lists of eigenvectors of graphs up to graph order, num steps for h
#has list of states sorted by order,
function generateListsofFidelityLists(order, graphs, hs, J, J2, width)
    list=Vector{Any}[];
    push!(list, getBaseFidelityList(order, graphs, hs));
    for i=1:order
        push!(list, generateFidelitiesList(order, graphs, hs, J, J2, width));
    end
    return list;
end

#here we dont do order+1 because we inputted the index inside the list, which is the correct index by itself
function getFidelityWeightAtParticularHValueandOrder(hIndex, order, fidelityLists, weights)
    sum=0;
    fid=fidelityLists[order][hIndex];
    for i=1:length(weights)
        sum+=weights[i];
    end
    #println("sum", sum)
    #println("fid", fid)
    #println("fid-sum", fid-sum);
    return fid-sum;
end

function getFidelityWeightsPerH(fidelityLists)
    #loop through all hs,
    weightsPerH=[];
    #fidelities list is a list of fidelities per order, each list having all the fidelities at different hs
    for hIndex=1:length(fidelityLists[1])
        weights=[];
        push!(weights, fidelityLists[1][hIndex])
        #loop through all the orders
        for order=2:length(fidelityLists)
            fidWeight=getFidelityWeightAtParticularHValueandOrder(hIndex, order, fidelityLists, weights)
            push!(weights, fidWeight);
        end
        push!(weightsPerH, weights);
    end
    return weightsPerH;
end

function calculateFidelitiesForEverything(order, graphs, J, J2, hs, width)
    hstep=(hs[length(hs)]-hs[1])/length(hs);
    statesList=getBaseStatesList(order, graphs, hs, J, J2, width, hstep);
    stateTemp=generateStatesList(order, graphs, hs, J, J2, width, hstep);

    append!(statesList, stateTemp);

    fidelities=[];

    for list in statesList
        fidTemp=calculateFidelity(hstep, list)
        push!(fidelities, fidTemp);
    end
    return fidelities;
end

function calculateInfiniteLatticeFidelity(order::Int, graphs, J, J2, hs, width)
    fidelitiesNLCResult=[];
    println("order, ", order);
    fidelities=calculateFidelitiesForEverything(order, graphs, J, J2, hs, width);
    println("fidelities ", fidelities);
    weights=getFidelityWeightsPerH(fidelities);
    println("weights ", weights);
    for i=1:length(weights)
        sum=weights[i][1];
        println("weights length", length(weights[i]));
        for j=2:length(weights[i])
            sum+=weights[i][j]*graphs[j-1].latticeConstant;
        end
        push!(fidelitiesNLCResult, sum);
    end
    return fidelitiesNLCResult;
end
#max order 56

function calculateInfiniteLatticeFidelity(orders::Vector{Int}, graphs, J, J2, hs, width)
    fidelitiesNLCResults=[];
    println("order ", orders[length(orders)]);
    fidelities=calculateFidelitiesForEverything(orders[length(orders)], graphs, J, J2, hs, width);
    println("fidelities, ", fidelities);
    weights=getFidelityWeightsPerH(fidelities);
    println("weights ", weights);
    for order in orders
        fid=[];
        for i=1:length(weights)
            println("length weights,should be order+1", length(weights[i]));
            println("order length", order+1);
            #println("order, ", order);
            sum=weights[i][1];
            for j=2:order+1
                println("length, ", length(weights[i]));
                sum+=weights[i][j]*graphs[j-1].latticeConstant;
            end
            push!(fid, sum);
        end
        push!(fidelitiesNLCResults, fid);
    end
    return fidelitiesNLCResults;
end


function calculateInfiniteLatticeEntanglement(orders::Vector{Int}, J, J2, h, graphs, width)
    @time begin
    println("weights starting!!");
    @time begin
        weights=getAllWeightsEntanglement(orders[length(orders)], graphs, J, J2, h, width);
    end

    println("weights done!!! ");
    ents=Any[];
    for order in orders
        sum=weights[1];
        for i=1:order
            #println("order: ", order);
            sum+=weights[i+1]*graphs[i].latticeConstant;
        end
        push!(ents, sum);
    end

end
    return ents;
end
