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
    println("eigenvector test, ", eigensystem[2]);
    sz=calculateSz(eigensystem[2], temp[3][eigensystem[3]], 1);
    return sz;
end

function getAllWeightsSz(num, graphs, J, J2, h, width)
    local weights::Vector{Float64}=Float64[];
    base=calculateBaseWeightSz(J, J2, h, width);
    push!(weights, base)
    for i=1:num
        #println("order: ", i);
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
    eigenvalues = temp[1];
    eigenvectors = temp[2];
    println("h: ", h);
    println("eigenvalues: ", eigenvalues);
    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    sz=theGraph.numSites*calculateSz(eigensystem[2], temp[3][eigensystem[3]], theGraph.numSites);
    for i=1:length(list)
        sum+=weights[list[i]+1];
    end
    #why???????? i thought weights[1] was the base weight with a site number of 1
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
    #why is this here??????
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
        bonds::Vector{bond}=bond[];
        temp=calculateEigensystemTransverseNoSymmetry(1, J, J2, h, bonds,"lanczos", "one", h, width, "H1");
        eigenvalues = temp[1]
        eigenvectors = temp[2]
        println("base state eigenvector: ", eigenvectors);
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



#=
function temp(order, graphs, J, J2, hs, width)
    hstep=(hs[length(hs)]-hs[1])/length(hs);
    statesList=getBaseStatesList(order, graphs, hs, J, J2, width, hstep);
    println(statesList);
    fidelity=calculateFidelity(hstep, statesList[1]);
    println("fidelity, ", fidelity);
end


=#

#SUSCEPTIBILITY NLC STARTS HERE
#N, J, J2, h, h2, bonds,eigmethod, num, hbar, hbar2, width, sites2
function calculateBaseWeightSusceptibility(J, J2, h, h2, width)
    bonds::Vector{bond}=bond[];
    #Eh calculation
    temp=calculateEigensystemSusceptibility(1, J, J2, h,h2, bonds,"lanczos", "one", h, h2, width, Int[1]);
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    #E0 calculation
    temp=calculateEigensystemTransverseNoSymmetry(1, J, J2, h, bonds,"lanczos", "one", h, width, "H1");
    eigenvalues1 = temp[1]
    eigenvectors1 = temp[2]
    Eh=eigenvalues[1];
    E0=eigenvalues1[1];
    sus=calculateSusceptibility(E0, Eh, h2);
    return sus;
end


function calculateWeightSusceptibility(num, graphs::Vector{graph}, weights, J, J2, h, h2, width)
    sum=0;
    #the graph to calculate weight of
    theGraph=graphs[num];
    list=theGraph.subgraphList;
    temp=copy(theGraph.nearBonds);
    bonds=append!(temp, theGraph.farBonds);
    #Eh calculation
    temp=calculateEigensystemSusceptibility(theGraph.numSites, J, J2, h, h2, bonds,"lanczos", "one", h, h2, width, theGraph.indicies);
    eigenvalues = temp[1];
    eigenvectors = temp[2];

    #E0 calculation
    temp=calculateEigensystemTransverse(theGraph.numSites, J, J2, h, bonds,"lanczos", "one", h, width);
    eigenvalues1 = temp[1];
    eigenvectors1 = temp[2];

    E0=eigenvalues1[1]
    Eh=eigenvalues[1]

    sus=theGraph.numSites*calculateSusceptibility(E0, Eh, h2);

    for i=1:length(list)
        sum+=weights[list[i]+1];
    end
    sum+=weights[1];
    return sus-sum;
end

function getAllWeightsSusceptibility(num, graphs, J, J2, h, h2, width)
    local weights::Vector{Float64}=Float64[];
    base=calculateBaseWeightSusceptibility(J, J2, h, h2, width);
    push!(weights, base)
    for i=1:num
        push!(weights, calculateWeightSusceptibility(i, graphs, weights, J, J2, h, h2, width));
    end
    return weights;
end

function calculateInfiniteLatticeSusceptibility(orders::Vector{Int}, J, J2, h, h2, graphs, width)
    @time begin
        sus=Any[];
        @time begin
            weights=getAllWeightsSusceptibility(orders[length(orders)], graphs, J, J2, h, h2, width);
        end

        for order in orders
            sum=weights[1];
            for i=1:order
                #println("order: ", order);
                sum+=weights[i+1]*graphs[i].latticeConstant;
            end
            push!(sus, sum);
        end

    end
    return sus;
end

#obtain a map mapping each site to the number of near bonds it has an number of far bonds
function obtainMapOfNumNearFarBonds(graph)
    list=Any[];
    nearBonds=graph.nearBonds;
    farBonds=graph.farBonds;
    near::Dict{Int,Int}=Dict{Int, Int}();
    far::Dict{Int,Int}=Dict{Int, Int}();
    for i=1:graph.numSites
        near[i]=0;
        far[i]=0;
    end
    for i=1:length(graph.farBonds)
        far[graph.farBonds[i].site1.num]+=1;
        far[graph.farBonds[i].site2.num]+=1;

    end
    for i=1:length(graph.nearBonds)
        near[graph.nearBonds[i].site1.num]+=1;
        near[graph.nearBonds[i].site2.num]+=1;

    end
    push!(list, near);
    push!(list, far);
    return list;
end


function obtainMapOfNumNearFarBonds(N, bonds)
    list=Any[];
    near::Dict{Int,Int}=Dict{Int, Int}();
    far::Dict{Int,Int}=Dict{Int, Int}();
    for i=1:N
        near[i]=0;
        far[i]=0;
    end

    for i=1:length(bonds)
        if(bonds[i].isNear)
            near[bonds[i].site1.num]+=1;
            near[bonds[i].site2.num]+=1;
        else
            far[bonds[i].site1.num]+=1;
            far[bonds[i].site2.num]+=1;
        end
    end
    push!(list, near);
    push!(list, far);
    return list;
end


#how many external sites is it connected to
#map1: how many far bonds connected to this site
#map2: how many near bonds connected
function obtainMeanFieldMapping(graph)
    #gives you a map of how many far bonds and near bonds each site has

    maps=obtainMapOfNumNearFarBonds(graph);
    map1=maps[1];
    map2=maps[2];

    near::Dict{Int,Int}=Dict{Int, Int}();
    far::Dict{Int,Int}=Dict{Int, Int}();
    list=Any[];

    for i=1:graph.numSites
        #random experiement
        factor=graph.indicies[i]==1 ? 1 : -1;
        near[i]=(4-map1[i])*factor;
        far[i]=-(2-map2[i])*factor;
    end
    push!(list, near);
    push!(list, far);
    return list;
end


function obtainMeanFieldMapping(N, bonds, maps, factors)
    map1=maps[1];
    map2=maps[2];

    near::Dict{Int,Int}=Dict{Int, Int}();
    far::Dict{Int,Int}=Dict{Int, Int}();
    list=Any[];

    for i=1:N
        near[i]=(4-map1[i])*factors[i];
        far[i]=-(2-map2[i])*factors[i];
    end
    push!(list, near);
    push!(list, far);
    return list;
end

function calculateActualField(site, map1, map2, h, J1, J2, m)
    #no more factors bc its integrated into map
    #maps are the factors
    value=map1[site]*m*J1+map2[site]*m*J2;
    return value;
end

#map 1: how many near bonds each site has
#map 2: how many far bonds each site has
#factors: whether the far/near bond has a factor of +-1.
function calculateSelfConsistentMz(numSites, map1, map2, h, J1, J2, bonds, J, firstGuess, maxIterations)
    hs=Float64[];
    ms=firstGuess;
    for i=1:numSites
        push!(hs, h);
    end
    if(ms==0)
        println("ERROR: don't input 0 initial guess");
        return 0;
    end
    count=0;
    percentError=1;
    while((abs(percentError)>0.00001)&&count<maxIterations)
        hs2=Float64[];
        println("iteration: ", count);
        for i=1:numSites
            push!(hs2,calculateActualField(i, map1, map2, h, J1, J2, ms));
        end
        temp=calculateEigensystemTransverseNoSymmetry(numSites, J1, J2, hs, hs2, bonds,"lanczos", "one", "H1");
        eigenvector = temp[2];
        states=temp[3][1];
        sz=calculateNeelOrderSz(eigenvector, states, numSites);
        percentError=(sz-ms);
        count+=1;
        println("initial guess: ", ms);
        println("calculated sz: ", sz);
        println("percentError", percentError);
        ms=sz;
    end
    return ms;
end

function calculateOneSelfConsistentSz(numSites, map1, map2, h, J1, J2, bonds, J, m, indicies)
    hs=Float64[];
    for i=1:numSites
        push!(hs, h);
    end
    if(m==0)
        println("ERROR: don't input 0 initial guess");
        return 0;
    end
    hs2=Float64[];
    for i=1:numSites
        push!(hs2,calculateActualField(i, map1, map2, h, J1, J2, m));
    end
    temp=calculateEigensystemTransverseNoSymmetry(numSites, J1, J2, hs, hs2, bonds,"lanczos", "one", "H1");
    eigenvector = temp[2];
    states=temp[3][1];
    sz=calculateNeelOrderSz(eigenvector, states, numSites, indicies);
    return sz;
end


function calculateSelfConsistentMz(numSites, map1, map2, h, J1, J2, bonds, J, firstGuess, maxIterations, indicies)
    hs=Float64[];
    ms=firstGuess;
    for i=1:numSites
        push!(hs, h);
    end
    if(ms==0)
        println("ERROR: don't input 0 initial guess");
        return 0;
    end
    count=0;
    percentError=1;
    #plot different iterations
    while((abs(percentError)>0.00001)&&count<maxIterations)
        hs2=Float64[];
        println("iteration: ", count);
        for i=1:numSites
            push!(hs2,calculateActualField(i, map1, map2, h, J1, J2, ms));
        end
        temp=calculateEigensystemTransverseNoSymmetry(numSites, J1, J2, hs, hs2, bonds,"lanczos", "one", "H1");
        eigenvector = temp[2];
        states=temp[3][1];
        sz=calculateNeelOrderSz(eigenvector, states, numSites, indicies);
        percentError=(sz-ms);
        count+=1;
        println("initial guess: ", ms);
        println("calculated sz: ", sz);
        println("percentError", percentError);
        ms=sz;
    end
    return ms;
end




function calculateBaseWeightMeanFieldSz(J, J2, h, firstGuess, maxIterations)
    bonds::Vector{bond}=bond[];
    map1::Dict{Int, Int}=Dict{Int, Int}();
    map2::Dict{Int, Int}=Dict{Int, Int}();
    #near bonds
    map1[1]=4;
    #far bonds
    map2[1]=2;
    sz=calculateOneSelfConsistentSz(1, map1, map2, h, J, J2, bonds, J, firstGuess, Int[1])
    println("Sz, ", sz);
    return sz;
end

#we are calculating spin in z direction so we want field in x direction, so H1
function calculateWeightMeanFieldSz(num, graphs::Vector{graph}, weights, J, J2, h, guess, maxIterations)
    sum=0;
    #the graph to calculate weight of
    theGraph=graphs[num];
    list=theGraph.subgraphList;
    maps=obtainMeanFieldMapping(theGraph);
    temp=copy(theGraph.nearBonds);
    bonds=append!(temp, theGraph.farBonds);
    sz=calculateOneSelfConsistentSz(theGraph.numSites, maps[1], maps[2], h, J, J2, bonds, J, guess, theGraph.indicies);
    for i=1:length(list)
        sum+=weights[list[i]+1];
    end
    println("this is sz for the ", num, "th graph: ", sz);
    sum+=theGraph.numSites*weights[1];
    return sz-sum;
end


function getAllWeightsMeanFieldSz(num, graphs, J, J2, h, firstGuess, maxIterations)
    local weights::Vector{Float64}=Float64[];
    #J, J2, h, firstGuess, maxIterations
    base=calculateBaseWeightMeanFieldSz(J, J2, h, firstGuess, maxIterations);
    push!(weights, base)
    for i=1:num
        #println("order: ", i);
        #num, graphs::Vector{graph}, weights, J, J2, h, firstGuess, maxIterations
        push!(weights, calculateWeightMeanFieldSz(i, graphs, weights, J, J2, h, firstGuess, maxIterations));
    end
    return weights;
end
#max order 56



function calculateInfiniteLatticeMeanFieldSz(orders::Vector{Int}, J, J2, h, graphs, firstGuess, maxIterations, absoluteError)
    @time begin
        #println("weights starting!!");
        szs=Any[];
        currGuess=Any[];
        #map is a map that maps a certain mean field value to a specific set of weights so that
        #you don't need to recalculate them
        map=Dict{Float64, Vector{Float64}}();
            weights=getAllWeightsMeanFieldSz(orders[length(orders)], graphs, J, J2, h, firstGuess, maxIterations);
            map[firstGuess]=weights;
        for i=1:length(orders)
            push!(currGuess, firstGuess);
        end

        #println("weights done!!! ");
        for order in orders
            sum=0;
            #println("order: ", order);
            ms=firstGuess;
            currWeights=map[ms];
            error=absoluteError+eps();
            count=0;
            while((abs(error)>absoluteError)&&count<maxIterations)
                sum=currWeights[1];
                for i=1:order
                    sum+=currWeights[i+1]*graphs[i].latticeConstant;
                end
                percentError=(sum-ms);
                count+=1;
                #println("initial guess: ", ms);
                #println("calculated sz: ", sum);
                #println("error", error);
                ms=sum;
                if(haskey(map, ms))
                    #wavefunction for a particular ms, should be the same????
                    currWeights= map[ms];
                else
                    map[ms]=getAllWeightsMeanFieldSz(orders[length(orders)], graphs, J, J2, h, ms, maxIterations);
                    currWeights=map[ms];
                    #println("the weights: ", currWeights);
                end
            end
            #so this is alllll the weights! so not everything was used... makes sense why there are 4...
            #must be something wrong with the second weight...
            println("the final weights: ", currWeights)
            push!(szs, sum);
        end
    end
    return szs;
end




function trackEachIterationSelfConsistentMz(numSites, map1, map2, h, J1, J2, bonds, J, firstGuess, maxIterations)
    hs=Float64[];
    ms=firstGuess;
    each=Any[];
    push!(each, ms);
    for i=1:numSites
        push!(hs, h);
    end
    if(ms==0)
        println("ERROR: don't input 0 initial guess");
        return 0;
    end
    count=0;
    while(count<maxIterations)
        hs2=Float64[];
        println("iteration: ", count);
        for i=1:numSites
            push!(hs2,calculateActualField(i, map1, map2, h, J1, J2, ms));
        end
        temp=calculateEigensystemTransverseNoSymmetry(numSites, J1, J2, hs, hs2, bonds,"lanczos", "one", "H1");
        eigenvector = temp[2];
        states=temp[3][1];
        sz=calculateNeelOrderSz(eigenvector, states, numSites);
        percentError=(sz-ms);
        count+=1;
        println("initial guess: ", ms);
        println("calculated sz: ", sz);
        println("percentError", percentError);
        ms=sz;
        push!(each, ms);
    end
    return each;
end
