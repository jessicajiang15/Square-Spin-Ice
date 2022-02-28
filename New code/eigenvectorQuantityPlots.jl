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

function plotSzVersusHTransverse(hList, J, J2, bonds, N)
    szlist=Float64[];
    for h=1:length(hList)
        temp = calculateEigensystemTransverse(N*N, J, J2, h, bonds,"lanczos", "one", h, 0);
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
        J2=0;
        hmin=0.2;
        hmax=1;
        num=50;
        deltah=(hmax-hmin)/num;
        hs=generateHListUniform(hmin, hmax, num);
        fids=Any[];
        bonds = bondListFrustrated(N)
        println("hs: ", hs);
        for i=1:length(hs)
            println("starting h: ", hs[i])
            temp = calculateEigensystemTransverse(N*N, J, J2, hs[i], bonds,"lanczos", "one", hs[i], 0);
            eigenvalues = temp[1]
            eigenvectors = temp[2]
            eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
            fid=calculateFidelity(eigensystem[2], temp[3][eigensystem[3]], hs[i], N, deltah, bonds, J, J2)
            #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
            push!(fids, fid);
        end
        println("fids: ", fids);
    end
    #TODO: plot it
    plot(hs, fids)
    savefig("./fidelityplot, J2 is 0.png")
end



function thesztest()
    println("Starting sz!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        J2=1;

        hs=generateHListUniform(0.1, 1, 20)
        ms=Any[];
        bonds = bondListFrustrated(N)
        println("hs: ", hs);
        for i=1:length(hs)
            println("starting h: ", hs[i])
            temp = calculateEigensystemTransverse(N*N, J, J2, hs[i], bonds,"lanczos", "one", hs[i], 0);

            eigenvalues = temp[1]
            eigenvectors = temp[2]
            eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);

            println("eigenvalues: ", eigenvalues);

            sz=calculateSz(eigensystem[2], temp[3][eigensystem[3]], N*N);
            #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
            push!(ms, sz);
            println("mz:", ms);
            println("sz:", sz);


        end
        println("ms: ", ms);

    end
    #TODO: plot it
    plot(hs, ms)
    savefig("./szplotnewnewnewlol ugh.png")
end



function calculateSpiTest()
    println("Starting sz!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        J2=0;
        hs=generateHListUniform(0.1, 2, 100);
        spis=Any[];
        bonds = bondListFrustrated(N)
        println("hs: ", hs);
        for i=1:length(hs)
            println("starting h: ", hs[i])
            temp =calculateEigensystemTransverseNoSymmetry(N*N, J, J2, hs[i], bonds,"lanczos", "one", hs[i], 0, "H1");
            eigenvalues = temp[1]
            eigenvectors = temp[2]
            #eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);

            #println("length of temp[3] ", length(temp[3]))
            #println("length of temp ", length(temp))
            #println("length of eigensystem ", length(eigensystem))

            #println("eigensystem[3] ", eigensystem[3]);
            println("starting time");
            @time begin
                spi=calculateSPiSzNew(eigenvectors[1], temp[3][1], N*N);
            end
            println("h: ", hs[i], " , spi: ", spi);

            push!(spis, spi);

        end
        println("spis: ", spis);
    end
    #TODO: plot it
    plot(hs, spis)
    savefig("./spiplot j2=0.png")
end

function calculateStaggeredFlippabilityTest()
    println("Starting sz!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        J2=1;

        hs=generateHListUniformHalf(J, J2, 50)
        flips=Any[];
        bonds = bondListFrustrated(N)
        println("hs: ", hs);
        squareIndicies=generateListsofPlaquetteIndiciesFlip(N);
        println("square", squareIndicies);
        for i=1:length(hs)
            println("starting h: ", hs[i])
            temp = calculateEigensystemTransverseNoSymmetry(N*N, J, J2, hs[i], bonds,"lanczos", "one", hs[i], 0, "H1");
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
        J2=1;

        order=20;
        hs=generateHListUniform(0.1, 0.5, 20)
        ms=Any[];
        graphs=readFromGraphFile();
        println("hs: ", hs);
        for i=1:length(hs)
            println("starting h: ", hs[i])
            sz=calculateInfiniteLatticeSz(order, J, J2, hs[i], graphs, 0)
            #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
            push!(ms, sz);
            println("sz:", sz);
        end
        println("ms: ", ms);

    end
    #TODO: plot it
    plot(hs, ms)
    savefig("./szplotinfinitelattice, order: "*string(order)*".png")
end






function weightsszgraphs()
    println("Starting sz inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        h=10;
        J2=1;

        order=2;
        width=0;
        ms=Any[];
        graphs=readFromGraphFile();
        weights=getAllWeightsSz(order, graphs, J, J2, h, width);
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
        J2=1;

        h=9.100000000000001;
        order=2;
        width=0;
        ms=Any[];
        graphs=readFromGraphFile();
        hs=generateHListUniform(1, 10, 10)
        println(hs);
        println("order, ", order)
        println("J, J2, ", J)
        println("h, ", h)
        println("width, ", width)
        weights=getAllWeightsSz(order, graphs, J, J2, hs[10], width);
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
        J2=1;

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
            println("J, J2, ", J)
            println("h, ", hs[i])
            println("width, ", width)
            weights=getAllWeightsSz(order, graphs, J, J2, hs[i], width);
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
    J2=1;

    width=0;
    ms=Any[];
    graphs=readFromGraphFile();
    hs=generateHListUniform(1, 10, 10)
    println(hs);
    j=length(hs)

    println("order, ", order)
    println("J, J2, ", J)
    println("width, ", width)
    println("h, ", hs[j])
    #println(graphs);
    #=
    weights=getAllWeightsSz(order, graphs, J, J2, hs[j-2], width);
    println("weights ", weights);
    =#
    weights1=getAllWeightsSz(order, graphs, J, J2, 2, width);
    #println("the fuck ", weights1);
    #println(graphs);
    #graphs=readFromGraphFile();

    weights=getAllWeightsSz(order, graphs, J, J2, hs[j], width);
    println("weights ", weights);


end




function weightsentanglementgraphs()
    println("Starting sz inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        J2=1;

        h=10;
        order=20;
        width=0;
        ms=Any[];
        graphs=readFromGraphFile();
        weights=getAllWeightsEntanglement(order, graphs, J, J2, h, width);
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
        J2=1;

        h=10;
        order=5;
        width=0;
        ms=Any[];
        graphs=readFromGraphFile();
        weights=getAllWeightsNoSubEntanglement(order, graphs, J, J2, h, width);
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
        J2=1;

        order=1;
        hs=generateHListUniform(0.1, 1, 50)
        ents=Any[];
        graphs=readFromGraphFile();
        println("hs: ", hs);
        for i=1:length(hs)
            println("starting h: ", hs[i])
            entanglement=calculateInfiniteLatticeEntanglement(order, J, J2, hs[i], graphs, 0)
            #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
            push!(ents, entanglement);
        end
        println("ents: ", ents);

    end
    #TODO: plot it
    plot(hs, ents)
    savefig("./entanglement NLC order: " * string(order)*".png")
end


function entanglementasj1j2()
    @time begin
        J2=0;
        J1=1;
        js=generateHListUniform(0.1, 2, 50);
        println(js);
        N=4;
        h=0.1;
        bonds = bondListFrustrated(N)
        listA=plaquetteIndicies(generateCheckerboardNoCrossPlaquettes(N)[1], N);
        println(listA);

        entropies=Any[];
        yes=false;
        for j in js
            if(j>1&&!yes)
                yes=true;
                temp = calculateEigensystemTransverse(N*N, J1, 1, h, bonds,"lanczos", "one", h, 0);
                eigenvalues = temp[1]
                eigenvectors = temp[2]
                eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
                entropy=getEntanglementEntropy(eigensystem[2], temp[3][eigensystem[3]], listA, N*N);
                push!(entropies, entropy);
                push!(js, 1);
            end
            println("j: ", j);
            temp = calculateEigensystemTransverse(N*N, J1, j, h, bonds,"lanczos", "one", h, 0);
            eigenvalues = temp[1]
            eigenvectors = temp[2]
            eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
            entropy=getEntanglementEntropy(eigensystem[2], temp[3][eigensystem[3]], listA, N*N);
            println(entropy);
            push!(entropies, entropy);
        end


        println(entropies);
        plot(js, entropies);
        savefig("./entanglemententropyplot"*", jmin: "*string(js[1])*", jmax: "*string(js[length(js)])*", h, "*string(h)*".png");
    end
end



function spij1j2()
    @time begin
        J1=1;
        js=generateHListUniform(0, 2, 400);
        println(js);
        N=4;
        h=0.1;
        bonds = bondListFrustrated(N)

        spis=Any[];

        for j in js
            println("j: ", j);
            temp =calculateEigensystemTransverseNoSymmetry(N*N, J1, j, h, bonds,"lanczos", "one", h, 0, "H1");
            eigenvalues = temp[1]
            eigenvectors = temp[2]
            spi=calculateSPiSzNew(eigenvectors[1], temp[3][1], N*N);
            println("spi, ", spi);
            push!(spis, spi);
        end
        println("spis, ", spis);
        plot(js, spis);
        savefig("./spiplot"*", jmin: "*string(js[1])*", jmax: "*string(js[length(js)])*", h: "*string(h)*".png");
    end
end





function weightspigraphs()
    println("Starting sz inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        J2=1;

        h=10;
        order=20;
        width=0;
        ms=Any[];
        graphs=readFromGraphFile();
        weights=getAllWeightsSpi(order, graphs, J, J2, h, width);
        println("weights", weights);

    end
    #TODO: plot it
    plot(1:length(weights), weights)
    savefig("./spisweightsinflattice.png")
end


function putin(list1, list2)
    for i=1:length(list2)
        push!(list1[i], list2[i]);
    end
end


function spiinfinitelattice()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        J2=0;
        os=Int[6];
        bonds = bondListFrustrated(N)
        hmin=0.1;
        hmax=2;
        hs=generateHListUniform(hmin, hmax, 100)
        println("hs: ", hs);
        graphs=readFromGraphFile();
        orders=Int[];
        list=Vector{Float64}[];
        spisother=Float64[];


        for i=1:length(os)
            push!(list, Vector{Float64}[]);
        end
        for i=1:length(os)
            push!(orders, getLastGraphNumOrder(os[i], graphs))
        end
        println("hs: ", hs);
        for i=1:length(hs)
            println("starting h: ", hs[i])
            @time begin
                spis=calculateInfiniteLatticeSpi(orders, J, J2, hs[i], graphs, 0)
                putin(list, spis);
                temp =calculateEigensystemTransverseNoSymmetry(N*N, J, J2, hs[i], bonds,"lanczos", "one", hs[i], 0, "H1");
                eigenvalues = temp[1]
                eigenvectors = temp[2]
                eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
                spi=calculateSPiSzNew(eigenvectors, temp[3][1], N*N);
                push!(spisother, spi);
            end
            #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
        end
    end
    #TODO: plot it
    push!(list, spisother);
    println("ED spis: ", spisother);
    plot(hs, list[1], label="order "*string(os[1]))

    println("spiinfinitelattice nlc: ", list);


    for i=2:length(list)
        str=i<=length(os) ? "order "*string(os[i]) : "ED";
        if(i!=length(list))
            plot!(hs, list[i], label=str)
        else
            plot!(hs, list[i], label=str)
        end
    end
    savefig("./spi NLC orders with ED, orders: " * string(os) *", j2: "*string(J2)*", hs: "*string(length(hs))*" hmin, "*string(hmin)*" hmax, "*string(hmax)*".png")
end


function multipleGraphsSz()
    maxOrder=5;


end


function spiinfinitelatticeone()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        J2=1;
        #hs=generateHListUniform(0.1, 1, 20)
        h=100
        orders=Int[];
        push!(orders, 22);
        graphs=readFromGraphFile();
        list=Vector{Float64}[];
        spi=calculateInfiniteLatticeSpi(orders, J, J2, h, graphs, 0)
        println("spi", spi);
    end
end




function spiinfinitelatticeold()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        J2=1;
        o=1;
        hs=generateHListUniform(0.1, 1, 20)
        spis=Any[];
        graphs=readFromGraphFile();
        order::Int=getLastGraphNumOrder(o, graphs);
        println("order: ", order);
        println("hs: ", hs);
        for i=1:length(hs)
            println("starting h: ", hs[i])
            @time begin
                spi=calculateInfiniteLatticeSpi(order, J, J2, hs[i], graphs, 0)
            end
            #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
            push!(spis, spi);
        end
        println("spis: ", spis);

    end
    #TODO: plot it
    plot(hs, spis)
    savefig("./spis NLC real order: " * string(o) *", "*".png")
end






function szinfinitelatticenew()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        J2=1;
        os=Int[1, 2, 3, 4, 5];
        hs=generateHListUniform(0.1, 1, 10)
        graphs=readFromGraphFile();

        orders=Int[];
        list=Vector{Float64}[];

        for i=1:length(os)
            push!(list, Vector{Float64}[]);
        end
        for i=1:length(os)
            push!(orders, getLastGraphNumOrder(os[i], graphs))
        end
        println("hs: ", hs);
        for i=1:length(hs)
            println("starting h: ", hs[i])
            @time begin
                szs=calculateInfiniteLatticeSz(orders, J, J2, hs[i], graphs, 0)
                putin(list, szs);
            end
            #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
        end
    end
    #TODO: plot it
    plot(hs, list[1], label="order "*string(os[1]))

    for i=2:length(os)
        plot!(hs, list[i], label="order "*string(os[i]))
    end
    savefig("./sz NEW PLS CHECK NLC orders: " * string(os) *", hs: "*string(length(hs))*".png")
end



function entanglementinfinitelatticenew()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        bonds = bondListFrustrated(N)
        J2=1;
        os=Int[4, 5, 6];
        #listA=plaquetteIndicies(generateCheckerboardPlaquettes(N)[1], N);
        listA=plaquetteIndicies(generateCheckerboardNoCrossPlaquettes(N)[1], N);

        hs=generateHListUniform(0.1, 1, 50)
        graphs=readFromGraphFile();
        entropies=Any[];

        orders=Int[];
        list=Vector{Float64}[];

        for i=1:length(os)
            push!(list, Vector{Float64}[]);
        end
        for i=1:length(os)
            push!(orders, getLastGraphNumOrder(os[i], graphs))
        end
        println("hs: ", hs);
        for i=1:length(hs)
            println("starting h: ", hs[i])
            @time begin
                ents=calculateInfiniteLatticeEntanglement(orders, J, J2, hs[i], graphs, 0)
                putin(list, ents);
                temp = calculateEigensystemTransverse(N*N, J, J2, hs[i], bonds,"lanczos", "one", hs[i], 0);
                eigenvalues = temp[1]
                eigenvectors = temp[2]
                eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
                entropy=getEntanglementEntropy(eigensystem[2], temp[3][eigensystem[3]], listA, N*N);
                push!(entropies, entropy);
            end
            #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
        end
    end
    #TODO: plot it
    push!(list, entropies);

    println("entropies ED: ", entropies);

    plot(hs, list[1].*2, label="order "*string(os[1]))

    println("entropies nlc: ", list);


    for i=2:length(list)
        str=i<=length(os) ? "order "*string(os[i]) : "ED";
        if(i!=length(list))
            plot!(hs, list[i].*2, label=str)
        else
            plot!(hs, list[i], label=str)
        end
    end
    savefig("./entanglement NLC orders: " * string(os) *", hs: "*string(length(hs))*".png")
end




function flippabilityj1j2()
    @time begin
        J1=1;
        js=generateHListUniform(0, 2, 100);
        println(js);
        N=4;
        h=0.1;
        bonds = bondListFrustrated(N)
        squareIndicies=generateListsofPlaquetteIndiciesFlip(N);

        flips=Any[];

        for j in js
            println("j: ", j);
            temp =calculateEigensystemTransverseNoSymmetry(N*N, J1, j, h, bonds,"lanczos", "one", h, 0, "H1");
            eigenvalues = temp[1]
            eigenvectors = temp[2]
            eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
            flip=calculateFlippabilityNew(eigenvectors[1], temp[3][1], N, squareIndicies)
            println("flip, ", flip);
            push!(flips, flip);
        end
        println("flips, ", flips);
        plot(js, flips);
        savefig("./flipplot"*", jmin: "*string(js[1])*", jmax: "*string(js[length(js)])*", h: "*string(h)*".png");
    end
end





function entanglementinfinitelatticenewj1j2()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        #listA=plaquetteIndicies(generateCheckerboardPlaquettes(N)[1], N);
        hs=generateHListUniform(0.1, 2, 100);
        js=generateHListUniformIncludeOne(0.1, 2, 5);
        graphs=readFromGraphFile();
        order=getLastGraphNumOrder(5, graphs);
        entropies=Any[];
        #list of different order 5s at different j values
        list=Vector{Float64}[];

        for j in js
            println("starting j: ", j)
            temp=Float64[];
            for i=1:length(hs)
                println("starting h: ", hs[i])
                @time begin
                    entanglement=calculateInfiniteLatticeEntanglement(order, J, j, hs[i], graphs, 0)
                    push!(temp, entanglement);
                end
                #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
            end
            push!(list, temp);
        end
    end
    println("entropies j1 j2: ", list);

    plot(hs, list[1], label="J2: "*string(js[1]))

    println("entropies nlc: ", list);


    for i=2:length(list)
        str="J2: "*string(js[i]);
        plot!(hs, list[i], label=str)
    end
    savefig("./entanglement NLC J2J1 vs h, order: " * string(order) *", hs: "*string(length(hs))*", js: "*string(length(js))*".png")
end




function entanglementEDj1j2vsh()
    println("Starting!!");
    @time begin
        N=4;
        J=1
        J2=0
        hmin=0.1;
        hmax=0.6;
        allData=Any[]
        hs=generateHListUniform(hmin, hmax, 50)
        println("hs", hs);
        js=generateHListUniformIncludeOne(0.9, 1.1, 10);
        listA=plaquetteIndicies(generateCheckerboardNoCrossPlaquettes(N)[1], N);
        println(listA);
        bonds = bondListFrustrated(N)
        for j in js
            entropies=Any[];
            for i=1:length(hs)
                temp = calculateEigensystemTransverse(N*N, J, j, hs[i], bonds,"lanczos", "one", hs[i], 0);
                eigenvalues = temp[1]
                eigenvectors = temp[2]
                eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
                entropy=getEntanglementEntropy(eigensystem[2], temp[3][eigensystem[3]], listA, N*N);
                push!(entropies, entropy);
            end
            push!(allData, entropies);
        end

        #println(entropy);
    end
    plot(hs, allData[1], label="j2: "*string(js[1]));
    for i=2:length(allData)
        plot!(hs, allData[i], label="j2: "*string(js[i]));
    end
    println(allData)
    savefig("./entanglementplot J2J1 js: "*string((js))*", hmin: "*string(hmin)*", hmax: "*string(hmax)*".png");

end




function calculateSpij1j2ED()
    println("Starting sz!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        hmin=0.1;
        hmax=2;
        J2=0;
        hs=generateHListUniform(hmin, hmax, 50);
        spis=Any[];
        bonds = bondListFrustrated(N)
        js=generateHListUniformIncludeOne(0.1, 2, 5);
        println("hs: ", hs);
        for j in js
            te=Any[];
            for i=1:length(hs)
                println("starting h: ", hs[i])
                temp =calculateEigensystemTransverseNoSymmetry(N*N, J, j, hs[i], bonds,"lanczos", "one", hs[i], 0, "H1");
                eigenvalues = temp[1]
                eigenvectors = temp[2]
                println("starting time");
                @time begin
                    spi=calculateSPiSzNew(eigenvectors[1], temp[3][1], N*N);
                end
                push!(te, spi);

            end
            push!(spis, te);

        end

        println("spis: ", spis);

        println("length spis", length(spis))
        println("length js", length(js));

        #println("spis: ", spis);
        plot(hs ,spis[1], label="j2: "*string(js[1]));
        for i=2:length(spis)
            plot!(hs ,spis[i], label="j2: "*string(js[i]));
        end
    end
    #TODO: plot it
    savefig("./spi J2J1 js: "*string((js))*", hmin: "*string(hmin)*", hmax: "*string(hmax)*".png");
end





function spiNLCj1j2()
    println("Starting entanglement inf!!");
    println("Starting sz!!")
    @time begin
        N=4;
        J=1
        J2=1;
        js=generateHListUniformIncludeOne(0.1, 2, 5);
        o=5;
        hs=generateHListUniform(0.1, 2, 100)
        spis=Any[];
        graphs=readFromGraphFile();
        order::Int=getLastGraphNumOrder(o, graphs);
        println("order: ", order);
        println("hs: ", hs);
        println("js: ", js);

        for j in js
            te=Any[];
            for i=1:length(hs)
                println("starting h: ", hs[i])
                @time begin
                    spi=calculateInfiniteLatticeSpi(order, J, j, hs[i], graphs, 0)
                end
                #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
                push!(te, spi);
            end
            push!(spis, te);
        end

        println("spis: ", spis);

    end
    #TODO: plot it
    plot(hs, spis[1], label="j2: "*string(js[1]));
    for i=2:length(spis)
        plot!(hs, spis[i], label="j2: "*string(js[i]));
    end
    savefig("./spis NLC j2j1 order: " * string(o) *", js: "*string(js)*", num h: "*string(length(hs))*".png")
end




function calculatefidelityJ1J2ED()
    println("Starting fid!!");
    @time begin
        N=4;
        J=1
        hmin=0.1;
        hmax=1;
        num=50;
        hs=generateHListUniform(hmin, hmax, num);
        js=generateHListUniformIncludeOne(0.1, 2, 5);
        fids=Any[];
        bonds = bondListFrustrated(N)
        println("hs: ", hs);
        for j in js
            te=calculateFidelity(hmin, hmax, num, N*N, J, j, bonds)
            push!(fids, te);
        end

        println("fids: ", fids);
    end
    #TODO: plot it
    plot(hs, fids[1], label="j2: "*string(js[1]))
    for i=2:length(fids)
        plot!(hs, fids[i], label="j2: "*string(js[i]));
    end
    savefig("./fidelity ED j2j1, js: "*string(js)*", num h: "*string(length(hs))*".png")
end



function thesztestj1j2ED()
    println("Starting sz!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        js=generateHListUniformIncludeOne(0, 2, 5);
        hs=generateHListUniform(0.1, 2, 100)
        ms=Vector{Float64}[];
        bonds = bondListFrustrated(N)
        println("hs: ", hs);
        for j in js
            te=Float64[];
            for i=1:length(hs)
                println("starting h: ", hs[i])
                temp = calculateEigensystemTransverse(N*N, J, j, hs[i], bonds,"lanczos", "one", hs[i], 0);
                eigenvalues = temp[1]
                eigenvectors = temp[2]
                eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
                println("eigenvalues: ", eigenvalues);
                sz=calculateSz(eigensystem[2], temp[3][eigensystem[3]], N*N);
                #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
                push!(te, sz);
            end
            push!(ms, te);
        end
        println("ms: ", ms);
    end
    #TODO: plot it
    plot(hs, ms[1], label="j2: "*string(js[1]))
    for i=2:length(ms)
        plot!(hs, ms[i], label="j2: "*string(js[i]));
    end
    savefig("./sz ED j2j1, js: "*string(js)*", num h: "*string(length(hs))*".png")
end





function szNLCj1j2()
    println("Starting sz inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        o=5;
        graphs=readFromGraphFile();

        order=getLastGraphNumOrder(o, graphs);
        hs=generateHListUniform(0.1, 2, 50)
        js=generateHListUniform(0.1, 2, 5)

        ms=Any[];
        println("hs: ", hs);
        for j in js
            te=Any[];
            for i=1:length(hs)
                println("starting h: ", hs[i])
                sz=calculateInfiniteLatticeSz(order, J, j, hs[i], graphs, 0)
                #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
                push!(te, sz);
            end
            push!(ms, te);
        end

        println("ms: ", ms);

    end
    #TODO: plot it
    plot(hs, ms[1], label="j2: "*string(js[1]))
    for i=2:length(ms)
        plot!(hs, ms[i], label="j2: "*string(js[i]));
    end
    savefig("./sz NLC j2j1, js: "*string(js)*", num h: "*string(length(hs))*", order: "*string(o)*".png")
end




function szinfinitelatticewithED()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        J2=1;
        os=Int[5, 6];
        hs=generateHListUniform(0.1, 1, 50)
        graphs=readFromGraphFile();
        bonds = bondListFrustrated(N)

        ms=Any[];

        orders=Int[];
        list=Vector{Float64}[];

        for i=1:length(os)
            push!(list, Vector{Float64}[]);
        end
        for i=1:length(os)
            push!(orders, getLastGraphNumOrder(os[i], graphs))
        end
        println("hs: ", hs);
        for i=1:length(hs)
            println("starting h: ", hs[i])
            @time begin
                szs=calculateInfiniteLatticeSz(orders, J, J2, hs[i], graphs, 0)
                putin(list, szs);
                temp = calculateEigensystemTransverse(N*N, J, J2, hs[i], bonds,"lanczos", "one", hs[i], 0);

                eigenvalues = temp[1]
                eigenvectors = temp[2]
                eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);

                println("eigenvalues: ", eigenvalues);

                sz=calculateSz(eigensystem[2], temp[3][eigensystem[3]], N*N);
                #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
                push!(ms, sz);
            end
            #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
        end
    end
    #TODO: plot it
    push!(list, ms);

    println(ms);

    plot(hs, list[1], label="order "*string(os[1]))

    for i=2:length(list)
        str=i<=length(os) ? "order "*string(os[i]) : "ED";
        if(i!=length(list))
            plot!(hs, list[i], label=str)
        else
            plot!(hs, list[i], label=str)
        end
    end
    savefig("./sz NLC with ED orders: " * string(os) *", hs: "*string(length(hs))*".png")
end


function spiinfinitelatticemultipleJsandorders()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        js=generateHListUniformIncludeOne(0.1, 1, 5);
        os=Int[1, 2, 3, 4, 5, 6];
        bonds = bondListFrustrated(N)
        hmin=9;
        hmax=10;
        hs=generateHListUniform(hmin, hmax, 2)
        println("hs: ", hs);
        graphs=readFromGraphFile();
        orders=Int[];


        for i=1:length(os)
            push!(orders, getLastGraphNumOrder(os[i], graphs))
        end
        println("hs: ", hs);
        for j in js
            spisother=Float64[];
            list=Vector{Float64}[];

            for i=1:length(os)
                push!(list, Vector{Float64}[]);
            end
            for i=1:length(hs)
                println("starting h: ", hs[i])
                spis=calculateInfiniteLatticeSpi(orders, J, j, hs[i], graphs, 0)
                putin(list, spis);
                temp =calculateEigensystemTransverseNoSymmetry(N*N, J, j, hs[i], bonds,"lanczos", "one", hs[i], 0, "H1");
                eigenvector = temp[2]

                spi=calculateSPiSzNew(eigenvector, temp[3][1], N*N);
                push!(spisother, spi);
                #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
            end
            #TODO: plot it
            push!(list, spisother);
            println("ED spis: ", spisother);
            plot(hs, list[1], label="order "*string(os[1]))

            println("spiinfinitelattice nlc: ", list);


            for i=2:length(list)
                str=i<=length(os) ? "order "*string(os[i]) : "ED";
                if(i!=length(list))
                    plot!(hs, list[i], label=str)
                else
                    plot!(hs, list[i], label=str)
                end
            end
            savefig("./spi NLC orders with ED, J2, " *string(j)*", orders, "* string(os) *", hs: "*string(length(hs))*" hmin, "*string(hmin)*" hmax, "*string(hmax)*".png")
        end
    end
end







function entanglementinfinitelatticemultipleJ2andorders()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        hmin=0.1;
        hmax=5;
        js=generateHListUniformIncludeOne(0.1, 2, 5)
        bonds = bondListFrustrated(N)
        os=Int[1, 2, 3, 4, 5];
        #listA=plaquetteIndicies(generateCheckerboardPlaquettes(N)[1], N);
        listA=plaquetteIndicies(generateCheckerboardNoCrossPlaquettes(N)[1], N);

        hs=generateHListUniform(hmin, hmax, 100)
        graphs=readFromGraphFile();

        orders=Int[];


        for i=1:length(os)
            push!(orders, getLastGraphNumOrder(os[i], graphs))
        end
        println("hs: ", hs);
        for j in js
            list=Vector{Float64}[];
            entropies=Any[];

            for i=1:length(os)
                push!(list, Vector{Float64}[]);
            end
            for i=1:length(hs)
                println("starting h: ", hs[i])
                @time begin
                    ents=calculateInfiniteLatticeEntanglement(orders, J, j, hs[i], graphs, 0)
                    putin(list, ents);
                    temp = calculateEigensystemTransverse(N*N, J, j, hs[i], bonds,"lanczos", "one", hs[i], 0);
                    eigenvalues = temp[1]
                    eigenvectors = temp[2]
                    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
                    entropy=getEntanglementEntropy(eigensystem[2], temp[3][eigensystem[3]], listA, N*N);
                    push!(entropies, entropy);
                end
            end
            push!(list, entropies);

            println("entropies ED: ", entropies);

            plot(hs, list[1].*2, label="order "*string(os[1]))

            println("entropies nlc: ", list);


            for i=2:length(list)
                str=i<=length(os) ? "order "*string(os[i]) : "ED";
                if(i!=length(list))
                    plot!(hs, list[i].*2, label=str)
                else
                    plot!(hs, list[i], label=str)
                end
            end
            savefig("./entanglement NLC orders with ED, J2, " *string(j)*", orders, "* string(os) *", hs: "*string(length(hs))*" hmin, "*string(hmin)*" hmax, "*string(hmax)*".png")
            #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
        end
    end
end







function szinfinitelatticemanyJ2andorders()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        js=generateHListUniformIncludeOne(0, 2, 5);
        os=Int[1, 2, 3, 4, 5];
        hmin=0.1;
        hmax=2.5;
        hs=generateHListUniform(hmin, hmax, 50)
        graphs=readFromGraphFile();
        bonds = bondListFrustrated(N)


        orders=Int[];

        for i=1:length(os)
            push!(orders, getLastGraphNumOrder(os[i], graphs))
        end
        println("hs: ", hs);

        for j in js
            list=Vector{Float64}[];
            ms=Any[];

            for i=1:length(os)
                push!(list, Vector{Float64}[]);
            end
            for i=1:length(hs)
                println("starting h: ", hs[i])
                @time begin
                    szs=calculateInfiniteLatticeSz(orders, J, j, hs[i], graphs, 0)
                    putin(list, szs);
                    temp = calculateEigensystemTransverse(N*N, J, j, hs[i], bonds,"lanczos", "one", hs[i], 0);

                    eigenvalues = temp[1]
                    eigenvectors = temp[2]
                    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
                    sz=calculateSz(eigensystem[2], temp[3][eigensystem[3]], N*N);
                    #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
                    push!(ms, sz);
                end
                #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
            end
            #TODO: plot it
            push!(list, ms);

            plot(hs, list[1], label="order "*string(os[1]))

            println(ms);

            for i=2:length(list)
                str=i<=length(os) ? "order "*string(os[i]) : "ED";
                if(i!=length(list))
                    plot!(hs, list[i], label=str)
                else
                    plot!(hs, list[i], label=str)
                end
            end
            savefig("./sz NLC orders with ED, J2, " *string(j)*", orders, "* string(os) *", hs: "*string(length(hs))*" hmin, "*string(hmin)*" hmax, "*string(hmax)*".png")
        end

    end

end


function fidelityInfiniteLattice()
    N=4;
    J=1
    hmin=0.1;
    hmax=2;
    num=20
    hs=generateHListUniform(hmin, hmax, num)
    println("num h: ", length(hs));
    graphs=readFromGraphFile();
    J2=0;
    width=0;
    bonds = bondListFrustrated(N)

    #te=calculateFidelity(hmin, hmax, num, N*N, J, J2, bonds);


    os=Int[1, 2, 3, 4];

    orders=Int[];

    for i=1:length(os)
        push!(orders, getLastGraphNumOrder(os[i], graphs))
    end

    println("orders ",orders);

    fidelities=calculateInfiniteLatticeFidelity(orders, graphs, J, J2, hs, width);
    println(fidelities);


    display(plot(hs[1:length(fidelities[1])], fidelities[1], title="NLC Fidelity susceptibility vs. h", label="Order "*string(os[1]), xlabel="h", ylabel="Fidelity"))

    println("fidelities nlc: ", fidelities);


    for i=2:length(fidelities)
        str="Order "*string(os[i]);
        display(plot!(hs[1:length(fidelities[i])], fidelities[i], label=str));
    end
    #display(plot!(hs, te, label="ED"))
    savefig("./fidelity NLC orders with ED, J2, " *string(J2)*", orders, "* string(os) *", hs: "*string(length(hs))*" hmin, "*string(hmin)*" hmax, "*string(hmax)*".png")
    #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
end




function fidelityInfiniteLattice2()
    N=4;
    J=1
    hmin=0.1;
    hmax=1;
    hs=generateHListUniform(hmin, hmax, 10)
    println("num h: ", length(hs));
    graphs=readFromGraphFile();

    J2=1;
    width=0;

    order=getLastGraphNumOrder(2, graphs);

    fidelities=calculateInfiniteLatticeFidelity(order, graphs, J, J2, hs, width);
    println("fidelities nlc: ", fidelities);
    te=calculateFidelity(hmin, hmax, num, N*N, J, j, bonds);


    display(plot(hs[1:length(fidelities)], fidelities, label="order "*string(order)))



end


function fidelityInfiniteLatticeMultipleJ2()
    N=4;
    J=1
    hmin=0.1;
    hmax=2;
    num=50;
    js=generateHListUniformIncludeOne(0.1, 1, 5);
    hs=generateHListUniform(hmin, hmax, num)
    println("num h: ", length(hs));
    graphs=readFromGraphFile();
    width=0;
    bonds = bondListFrustrated(N)

    os=Int[1, 2, 3, 4, 5];

    orders=Int[];

    for i=1:length(os)
        push!(orders, getLastGraphNumOrder(os[i], graphs))
    end

    println("orders ",orders);
    for j in js
        fidelities=calculateInfiniteLatticeFidelity(orders, graphs, J, j, hs, width);
        println("j2=", j, ", fidelities: ", fidelities);
        display(plot(hs[1:length(fidelities[1])], fidelities[1], title="NLC Fidelity susceptibility vs. h, J2="*string(j), label="Order "*string(os[1]), xlabel="h", ylabel="Fidelity"))
        println("fidelities nlc: ", fidelities);
        te=calculateFidelity(hmin, hmax, num, N*N, J, j, bonds);
        for i=2:length(fidelities)
            str="Order "*string(os[i]);
            display(plot!(hs[1:length(fidelities[i])], fidelities[i], label=str));
        end
        display(plot!(hs, te, label="ED"))
        savefig("./fidelity NLC orders with ED new, J2, " *string(j)*", orders, "* string(os) *", hs: "*string(length(hs))*" hmin, "*string(hmin)*" hmax, "*string(hmax)*".png")
    end
end

function susceptibilityMultipleJ2J1AndOrders()
    println("Starting sus inf!!");
    @time begin
        J=1
        N=4;
        js=generateHListUniformIncludeOne(0.1, 1, 5);
        os=Int[1, 2, 3, 4, 5, 6];
        bonds = bondListFrustrated(N)
        hmin=9;
        h2=0.1;
        hmax=10;
        hs=generateHListUniform(hmin, hmax, 2)
        println("hs: ", hs);
        graphs=readFromGraphFile();
        orders=Int[];
        nums=getPlaquetteNumberList(4);


        for i=1:length(os)
            push!(orders, getLastGraphNumOrder(os[i], graphs))
        end
        println("hs: ", hs);
        for j in js
            suses=Float64[];
            list=Vector{Float64}[];

            for i=1:length(os)
                push!(list, Vector{Float64}[]);
            end
            for i=1:length(hs)
                println("starting h: ", hs[i])
                sus=calculateInfiniteLatticeSusceptibility(orders, J, j, hs[i], h2, graphs, 0)
                putin(list, sus);
                temp =calculateEigensystemTransverseNoSymmetry(N*N, J, j, hs[i], bonds,"lanczos", "one", hs[i], 0, "H1");
                E0 = temp[1][1];
                temp2=calculateEigensystemSusceptibility(N, J, J2, h, h2, bonds,"lanczos", "one", h, h2, 0, nums);
                Eh=temp[1][1];
                su=calculateSusceptibility(E0, Eh, h2);
                push!(suses, su);
                #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
            end
            #TODO: plot it
            push!(list, suses);
            plot(hs, list[1], label="Order "*string(os[1]))

            println("susceptibilities nlc: ", list);


            for i=2:length(list)
                str=i<=length(os) ? "Order "*string(os[i]) : "Exact Diagonalization";
                if(i!=length(list))
                    plot!(hs, list[i], label=str)
                else
                    plot!(hs, list[i], label=str)
                end
            end
            savefig("./sus NLC orders with ED, J2, " *string(j)*", orders, "* string(os) *", hs: "*string(length(hs))*" hmin, "*string(hmin)*" hmax, "*string(hmax)*".png")
        end
    end
end






function susceptibilityJ2J1ED()
    println("Starting sz!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        nums=getPlaquetteNumberList(N);
        h2=0.1;
        js=generateHListUniformIncludeOne(0.1, 2, 5);
        println("js: ", js);
        hs=generateHListUniform(0.1, 2, 20)
        ms=Vector{Float64}[];
        bonds = bondListFrustrated(N);
        println("hs: ", hs);
        for j in js
            te=Float64[];
            for i=1:length(hs)
                println("starting h: ", hs[i])
                temp =calculateEigensystemTransverseNoSymmetry(N*N, J, j, hs[i], bonds,"lanczos", "one", hs[i], 0, "H1");
                E0 = temp[1][1];
                temp2=calculateEigensystemSusceptibility(N*N, J, j, hs[i], h2, bonds,"lanczos", "one", hs[i], h2, 0, nums);
                Eh=temp2[1][1];
                su=calculateSusceptibility(E0, Eh, h2);
                #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
                push!(te, su);
            end
            push!(ms, te);
        end
        println("sus: ", ms);
    end
    #TODO: plot it
    plot(hs, ms[1], label="J2: "*string(js[1]))
    for i=2:length(ms)
        plot!(hs, ms[i], label="J2: "*string(js[i]), xlabel="h", ylabel="Susceptibility");
    end
    savefig("./susceptibility ED js: "*string(js)*", num h: "*string(length(hs))*".png")
end




function susinfinitelatticemanyJ2andorders()
    println("Starting sus inf!!");
    println("Starting sus!!");

    @time begin
        N=4;
        h2=0.1;
        J=1
        js=generateHListUniformIncludeOne(0.1, 2, 5);
        os=Int[1, 2, 3, 4, 5];
        hmin=0.1;
        nums=getPlaquetteNumberList(N);

        hmax=2;
        hs=generateHListUniform(hmin, hmax, 50)
        graphs=readFromGraphFile();
        bonds = bondListFrustrated(N)


        orders=Int[];

        for i=1:length(os)
            push!(orders, getLastGraphNumOrder(os[i], graphs))
        end
        println("hs: ", hs);

        for j in js
            list=Vector{Float64}[];
            ms=Any[];

            for i=1:length(os)
                push!(list, Vector{Float64}[]);
            end
            for i=1:length(hs)
                println("starting h: ", hs[i])
                @time begin
                    szs=calculateInfiniteLatticeSusceptibility(orders, J, j, hs[i], h2, graphs, 0)
                    putin(list, szs);
                    temp =calculateEigensystemTransverseNoSymmetry(N*N, J, j, hs[i], bonds,"lanczos", "one", hs[i], 0, "H1");
                    E0 = temp[1][1];
                    temp2=calculateEigensystemSusceptibility(N*N, J, j, hs[i], h2, bonds,"lanczos", "one", hs[i], h2, 0, nums);
                    Eh=temp2[1][1];
                    sz=calculateSusceptibility(E0, Eh, h2);
                    #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
                    push!(ms, sz);
                end
                #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
            end
            #TODO: plot it
            push!(list, ms);

            plot(hs, list[1], label="order "*string(os[1]))

            println(ms);

            for i=2:length(list)
                str=i<=length(os) ? "order "*string(os[i]) : "ED";
                if(i!=length(list))
                    plot!(hs, list[i], label=str)
                else
                    plot!(hs, list[i], label=str)
                end
            end
            savefig("./sus NLC orders with ED, J2, " *string(j)*", orders, "* string(os) *", hs: "*string(length(hs))*" hmin, "*string(hmin)*" hmax, "*string(hmax)*".png")
        end

    end

end





function szMeanFieldInfiniteLattice()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        J2=0;
        os=Int[1, 2, 3];
        hs=generateHListUniform(0.1, 1, 50)
        graphs=readFromGraphFile();
        firstGuess=0.1;
        maxIterations=30;

        orders=Int[];
        list=Vector{Float64}[];

        for i=1:length(os)
            push!(list, Vector{Float64}[]);
        end
        for i=1:length(os)
            push!(orders, getLastGraphNumOrder(os[i], graphs))
        end
        println("hs: ", hs);
        for i=1:length(hs)
            println("starting h: ", hs[i])
            @time begin
                #calculateInfiniteLatticeMeanFieldSz(orders::Vector{Int}, J, J2, h, graphs, firstGuess, maxIterations)
                szs=calculateInfiniteLatticeMeanFieldSz(orders, J, J2, hs[i], graphs, firstGuess, maxIterations)
                putin(list, szs);
            end
            #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
        end
    end
    #TODO: plot it
    plot(hs, list[1], label="order "*string(os[1]))

    for i=2:length(os)
        plot!(hs, list[i], label="order "*string(os[i]))
    end
    savefig("./sz MEAN FIELD orders: " * string(os) *", hs: "*string(length(hs))*".png")
end






function meanFieldSzInfiniteLatticeManyJ2AndOrders()
    println("Starting entanglement inf!!");
    println("Starting sz!!");

    @time begin
        N=4;
        J=1
        js=[0]
        #js=generateHListUniform(0.1, 1, 5);
        os=Int[1,2];
        hmin=0.1;
        hmax=0.11;
        hs=generateHListUniform(hmin, hmax, 2)
        graphs=readFromGraphFile();
        bonds = bondListFrustrated(N)
        firstGuess=0.1;
        maxIterations=20;
        absError=0.0001;
        orders=Int[];
        for i=1:length(os)
            push!(orders, getLastGraphNumOrder(os[i], graphs))
        end
        println("hs: ", hs);
        println("orders: ", orders);

        for j in js
            list=Vector{Float64}[];
            ms=Any[];

            for i=1:length(os)
                push!(list, Vector{Float64}[]);
            end
            for i=1:length(hs)
                println("starting h: ", hs[i])
                @time begin
                    szs=calculateInfiniteLatticeMeanFieldSz(orders, J, j, hs[i], graphs, firstGuess, maxIterations, absError)
                    putin(list, szs);
                end
                #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
            end
            #TODO: plot it
            plot(hs, list[1], label="order "*string(os[1]))

            println(list);

            for i=2:length(list)
                str="order "*string(os[i]);
                plot!(hs, list[i], label=str)
            end
            savefig("./sz mean field calc NLC orders with ED THIS ONE, J2, " *string(j)*", orders, "* string(os) *", hs: "*string(length(hs))*" hmin, "*string(hmin)*" hmax, "*string(hmax)*".png")
        end

    end

end


function calculateEDMeanFieldSz()
    N=4;
    bonds = bondListFrustratedNoPBC(N)
    maps=obtainMapOfNumNearFarBonds(N*N, bonds);
    factors=calculateExternalFieldFactors(N, maps[1], maps[2]);
    result=obtainMeanFieldMapping(N*N, bonds, maps, factors);
    map1=result[1];
    map2=result[2];
    firstGuess=0.1;
    maxIterations=60;
        @time begin
            J=1
            js=[0];
            hs=generateHListUniform(0.1, 2, 50)
            ms=Vector{Float64}[];
            println("hs: ", hs);
            for j in js
                te=Float64[];
                for i=1:length(hs)
                    println("starting h: ", hs[i])
                    sz=calculateSelfConsistentMz(N*N, map1, map2, hs[i], J, j, bonds, J, firstGuess, maxIterations)
                    push!(te, sz);
                end
                push!(ms, te);
            end
            println("ms: ", ms);
        end
        #TODO: plot it
        plot(hs, ms[1], label="j2: "*string(js[1]))
        for i=2:length(ms)
            plot!(hs, ms[i], label="j2: "*string(js[i]));
        end
        savefig("./sz mean field ED j2j1, js: "*string(js)*", num h: "*string(length(hs))*".png")
end





function trackEDMeanFieldSzEachIteration()
    N=4;
    bonds = bondListFrustratedNoPBC(N)
    maps=obtainMapOfNumNearFarBonds(N*N, bonds);
    factors=calculateExternalFieldFactors(N, maps[1], maps[2]);
    result=obtainMeanFieldMapping(N*N, bonds, maps, factors);
    map1=result[1];
    map2=result[2];
    firstGuess=0.1;
    maxIterations=100;
        @time begin
            J=1
            J2=0.64;
            hmin=0.8
            hmax=1
            hs=generateHListUniform(hmin, hmax, 20);
            ms=Vector{Float64}[];
            println("hs: ", hs);
            te=Vector{Any}[];
                for i=1:length(hs)
                    println("starting h: ", hs[i])
                    each=trackEachIterationSelfConsistentMz(N*N, map1, map2, hs[i], J, J2, bonds, J, firstGuess, maxIterations)
                    push!(te, each);
                end
            println("ms: ", te);
        end
        #TODO: plot it
        plot(1:maxIterations+1, te[1], label="h: "*string(hs[1]))
        for i=2:length(te)
            plot!(1:maxIterations+1, te[i], label="h: "*string(hs[i]));
        end
        savefig("./sz mean field iterations ED, js: "*string(J2)*", num h: "*string(length(hs))*"hmin, "*string(hmin)*", hmax, "*string(hmax)*".png")
end
