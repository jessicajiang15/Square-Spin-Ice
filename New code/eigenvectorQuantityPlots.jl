include("hamiltonianDiagonalization.jl")
include("calculations and plotting.jl")



function plotSxVersusH(hList)

end


function findLowest(list)
    min=1;
    for i=2:length(list)
        if(list[i]<list[min])
            min=i;
        end
    end
    return min;
end

function plotSzVersusHTransverse(hList, J, bonds, N)
    szlist=Float64[];
    for h=1:length(hList)
        temp = calculateEigensystemTransverse(N*N, J, h, bonds,"lanczos", "one", h, 0);
        minInd=findLowest(temp[1]);
        eigenvector = temp[2][minInd];
        #temp[3] is a map mapping index to the list of states
        states= minInd>=temp[3][1][1]&&minInd<=temp[3][1][2] ? temp[3][2] : temp[3][4];
        push!(szlist, calculateSz(eigenvector, states, N));
    end
    return szlist;
end




function calculatefidelity()
    println("Starting fid!!");
    @time begin
    N=4;
    J=1
    deltah=0.001;
    hs=generateHListUniform(0.1, 1, 50);
    fids=Any[];
    bonds = bondListFrustrated(N)
    println("hs: ", hs);
    for i=1:length(hs)
        println("starting h: ", hs[i])
        temp = calculateEigensystemTransverse(N*N, J, hs[i], bonds,"lanczos", "one", hs[i], 0);
        eigenvalues = temp[1]
        eigenvectors = temp[2]
        eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
        fid=calculateFidelity(eigensystem[2], temp[3][eigensystem[3]], hs[i], N, deltah, bonds, J)
        #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
        push!(fids, fid);
    end
    println("fids: ", fids);
    end
    #TODO: plot it
    plot(hs, fids)
    savefig("./fidelityplot.png")

end



function thesztest()
    println("Starting sz!!");
    println("Starting sz!!");

    @time begin
    N=4;
    J=1
    hs=generateHListUniformHalf(J, 50)
    ms=Any[];
    bonds = bondListFrustrated(N)
    println("hs: ", hs);
    for i=1:length(hs)
        println("starting h: ", hs[i])
        temp = calculateEigensystemTransverse(N*N, J, hs[i], bonds,"lanczos", "one", hs[i], 0);
        eigenvalues = temp[1]
        eigenvectors = temp[2]
        eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);

        println("length of temp[3] ", length(temp[3]))
        println("length of temp ", length(temp))
        println("length of eigensystem ", length(eigensystem))
        println("eigensystem[3] ", eigensystem[3]);
        sz=calculateSz(eigensystem[2], temp[3][eigensystem[3]], N);
        #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
        push!(ms, sz);
        println("mz:", ms);
        println("sz:", sz);


    end
    println("ms: ", ms);

    end
    #TODO: plot it
    plot(hs, ms)
    savefig("./szplotnewnewnewlol.png")
end



function calculateSpiTest()
    println("Starting sz!!");
    println("Starting sz!!");

    @time begin
    N=4;
    J=1
    hs=generateHListUniform(0.1, 1, 100);
    spis=Any[];
    bonds = bondListFrustrated(N)
    println("hs: ", hs);
    for i=1:length(hs)
        println("starting h: ", hs[i])
        temp =calculateEigensystemTransverseNoSymmetry(N*N, J, hs[i], bonds,"lanczos", "one", hs[i], 0, "H1");
        eigenvalues = temp[1]
        eigenvectors = temp[2]
        #eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);

        #println("length of temp[3] ", length(temp[3]))
        #println("length of temp ", length(temp))
        #println("length of eigensystem ", length(eigensystem))

        #println("eigensystem[3] ", eigensystem[3]);
        println("starting time");
        @time begin
        spi=calculateSPiSzNew(eigenvectors[1], temp[3][1], N);
    end
        println("h: ", hs[i], " , spi: ", spi);

        push!(spis, spi);

    end
    println("spis: ", spis);
    end
    #TODO: plot it
    plot(hs, spis)
    savefig("./spiplotsmallh.png")
end

function calculateStaggeredFlippabilityTest()
    println("Starting sz!!");
    println("Starting sz!!");

    @time begin
    N=4;
    J=1
    hs=generateHListUniformHalf(J, 50)
    flips=Any[];
    bonds = bondListFrustrated(N)
    println("hs: ", hs);
    squareIndicies=generateListsofPlaquetteIndiciesFlip(N);
    println("square", squareIndicies);
    for i=1:length(hs)
        println("starting h: ", hs[i])
        temp = calculateEigensystemTransverseNoSymmetry(N*N, J, hs[i], bonds,"lanczos", "one", hs[i], 0, "H1");
        eigenvalues = temp[1]
        eigenvectors = temp[2]
        eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);

        flip=calculateFlippabilityNew(eigenvectors[1], temp[3][1], N, squareIndicies)
        println("h: ", hs[i], " , flip: ", flip);

        push!(flips, flip);
    end
    println("flips: ", flips);
    end
    #TODO: plot it
    plot(hs, flips)
    savefig("./flippabilitytest1.png")
end



function thesztestinfinitelattice()
    println("Starting sz inf!!");
    println("Starting sz!!");

    @time begin
    N=4;
    J=1
    order=20;
    hs=generateHListUniform(0.1, 1, 50)
    ms=Any[];
    graphs=readFromGraphFile();
    println("hs: ", hs);
    for i=1:length(hs)
        println("starting h: ", hs[i])
        sz=calculateInfiniteLatticeSz(order, J, hs[i], graphs, 0)
        #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
        push!(ms, sz);
        println("sz:", sz);
    end
    println("ms: ", ms);

    end
    #TODO: plot it
    plot(hs, ms)
    savefig("./szplotinfinitelatticesmallh.png")
end






function weightsszgraphs()
    println("Starting sz inf!!");
    println("Starting sz!!");

    @time begin
    N=4;
    J=1
    h=10;
    order=2;
    width=0;
    ms=Any[];
    graphs=readFromGraphFile();
        weights=getAllWeightsSz(order, graphs, J, h, width);
        println(weights);

    end
    #TODO: plot it
    plot(1:length(weights), weights)
    savefig("./szinfiniteweightsloworder.png")
end



function weightsszgraphs2()
    println("Starting sz inf!!");
    println("Starting sz!!");

    @time begin
    N=4;
    J=1
    h=9.100000000000001;
    order=2;
    width=0;
    ms=Any[];
    graphs=readFromGraphFile();
    hs=generateHListUniform(1, 10, 10)
    println(hs);
    println("order, ", order)
    println("J, ", J)
    println("h, ", h)
    println("width, ", width)
        weights=getAllWeightsSz(order, graphs, J, hs[10], width);
        println("weights ", weights);
    end
    #TODO: plot it
    println(weights);
    plot(1:length(weights), weights)
    savefig("./graph0weights1.png")


end



function weightsszgraphs3()
    println("Starting sz inf!!");
    println("Starting sz!!");

    @time begin
    J=1
    order=2;
    width=0;
    ms=Any[];
    graphs=readFromGraphFile();
    hs=generateHListUniform(1, 10, 10)
    list1=Float64[];
    list2=Float64[];
    list3=Float64[];
    for i=1:length(hs)
        println("order, ", order)
        println("J, ", J)
        println("h, ", hs[i])
        println("width, ", width)
        weights=getAllWeightsSz(order, graphs, J, hs[i], width);
        println("weights", weights);
        push!(list1, weights[1])
        push!(list2, weights[2])
        push!(list3, weights[3])
    end

    println("list2", list2)
    println("list3", list3);

    end
    #TODO: plot it
    plot(hs, list1)
    savefig("./graph0weights1.png")
    plot(hs, list2)
    savefig("./graph1weights1.png")

    plot(hs, list3)
    savefig("./graph2weights1.png")

end








function weightsszgraphs4()
    println("Starting sz inf!!");
    println("Starting sz!!");

    N=4;
    J=1
    order=2;
    width=0;
    ms=Any[];
    graphs=readFromGraphFile();
    hs=generateHListUniform(1, 10, 10)
    println(hs);
    j=length(hs)

        println("order, ", order)
        println("J, ", J)
        println("width, ", width)
        println("h, ", hs[j])
        #println(graphs);
#=
        weights=getAllWeightsSz(order, graphs, J, hs[j-2], width);
        println("weights ", weights);
=#
        weights1=getAllWeightsSz(order, graphs, J, 2, width);
        #println("the fuck ", weights1);
        #println(graphs);
        #graphs=readFromGraphFile();

        weights=getAllWeightsSz(order, graphs, J, hs[j], width);
        println("weights ", weights);


end




function weightsentanglementgraphs()
    println("Starting sz inf!!");
    println("Starting sz!!");

    @time begin
    N=4;
    J=1
    h=10;
    order=5;
    width=0;
    ms=Any[];
    graphs=readFromGraphFile();
        weights=getAllWeightsEntanglement(order, graphs, J, h, width);
        println("weights", weights);

    end
    #TODO: plot it
    plot(1:length(weights), weights)
    savefig("./entanglementinfiniteweightsloworder.png")
end


function weightsentanglementgraphsnosub()
    println("Starting sz inf!!");
    println("Starting sz!!");

    @time begin
    N=4;
    J=1
    h=10;
    order=5;
    width=0;
    ms=Any[];
    graphs=readFromGraphFile();
        weights=getAllWeightsNoSubEntanglement(order, graphs, J, h, width);
        println("weights", weights);
    end
    #TODO: plot it
    plot(1:length(weights), weights)
    savefig("./entanglementinfiniteweightsnosub.png")
end




function theentanglementtestinfinitelattice()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
    N=4;
    J=1
    order=5;
    hs=generateHListUniform(0.1, 10, 50)
    ents=Any[];
    graphs=readFromGraphFile();
    println("hs: ", hs);
    for i=1:length(hs)
        println("starting h: ", hs[i])
        entanglement=calculateInfiniteLatticeEntanglement(order, J, hs[i], graphs, 0)
        #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
        push!(ents, entanglement);
    end
    println("ents: ", ents);

    end
    #TODO: plot it
    plot(hs, ents)
    savefig("./entanglementplotinfinitelattice.png")
end
