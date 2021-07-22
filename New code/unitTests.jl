include("entanglement helper.jl")


function test()
    println("lmaoas");
    N=4;
    J=1
    h=1;
    hs=generateHListUniform(J, 20);
    ms=Any[];
    bonds = bondListFrustrated(N)
    temp = calculateEigensystemTransverse(N, J, h, bonds,"lanczos", "one", h, 0);
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
    sz=calculateSz(eigensystem[2], temp[3][eigensystem[3]], N);
    println("sz: ", sz);
end

function testMeasures()
    N=4;
    println("hello");
    states::Vector{Int}=[1, 2, 3];
    map::Dict{Int, Int}=Dict([(1, 1), (2, 2), (3, 3)]);
    testState::Vector{Float64}=[1/sqrt(2), 1/sqrt(4), 1/sqrt(4)];
    sz=calculateSz(testState, states, N);
    sx=calculateSx(testState, states, map, N);
    ind=indiciesAllSquares(N);
    stflippability=calculateStaggeredFlippability(testState, states, ind, N);
    println("sz, ", sz);
    println("sx, ", sx);
    println("flip, ", stflippability);

end

function testPlot()
    J=1;
    N=4;
    hList=generateHListUniform(J, 20);
    bonds = bondListFourNeighbors(N);
    plotSzVersusHTransverse(hList, J, bonds, N);
end

function testbond()
    bonds=bondListFrustrated(4)
    for i in bonds
        if(containsSite(1, i))
            println(i);
        end
    end
end

function thetest()
    #list1=Int[1, 2, 3];
    #list2=Int[1, 2, 3];
    x = 1:10; y = rand(10, 2) # 2 columns means two lines
    plot(x, y)

    list=Any[x, y];

    #println(typeof(list[1]));
    #theplot(list, 2, "lmao")

end

function entanglementtest()
    println("starting");
    println("starting");

    #00,01,10,11
    #0000, 0001, 0000.... 1111
    N=2
    eigenvector=Float64[1/sqrt(2), 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/sqrt(2)];
    println(length(eigenvector));
    listA=Int[0, 1];
    states::Array{Int}=Int[];
    for i=0:2^(N*N)-1
        push!(states, i);
    end
    println("calculating entropy");
    #eigenvector, states, listA, N
    println(getEntanglementEntropy(eigenvector, states, listA, N));
end

function lol()
        println("starting");
    println("lmao");
    N=4;
    listA=Int[0, 1];
    #=
    N=2;
    listA=Int[0, 1];
    for i=0:15
        println(getBNum(i, listA, 2))
    end
=#
    states::Array{Int}=Int[];
    for i=0:2^(N*N)-1
        push!(states, i);
    end

    for i=1:length(states)
        println("i, ",states[i]);
        getA=getANum(states[i], listA);
        getB=getBNum(states[i], listA, 2);
    end


end

function testExtractDigits()
    listA=Int[1, 2, 4, 6, 9, 10, 11, 12];
    getA=getANum(7766, listA);
    getB=getBNum(7766, listA, 2);
    println("A ", getA);
    println("B ", getB);
end

function calculategsentanglement()
    println("Starting!!");
    @time begin
    N=4;
    J=1
    entropies=Any[];

        hs=generateHListUniform(J, 20);
        listA=plaquetteIndicies(generateCheckerboardPlaquettes(N)[1], N);
        bonds = bondListFrustrated(N)
        for i=1:length(hs)
            temp = calculateEigensystemTransverse(N, J, hs[i], bonds,"lanczos", "one", hs[i], 0);
            eigenvalues = temp[1]
            eigenvectors = temp[2]
            eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
            entropy=getEntanglementEntropy(eigensystem[2], temp[3][eigensystem[3]], listA, N);
            push!(entropies, entropy);
        end
    #println(entropy);
    end
    println(entropies)
    plot(hs, entropies)
    #savefig("Users/Jessica/git/Square Spin Ice/szplot.png")
    savefig("./entanglementplot.png")

end


function graphTest()
x = [1, 2, 3]
y=[1,2,3]
plot(x, y)
display(plot(x, y))
savefig("./plot.png")
end


function thesztest()
    println("Starting sz!!");
    @time begin
    N=4;
    J=1
    hs=generateHListUniform(J, 20);
    ms=Any[];
    bonds = bondListFrustrated(N)
    println("hs: ", hs);
    for i=1:length(hs)
        println("starting h: ", hs[i])
        temp = calculateEigensystemTransverse(N, J, hs[i], bonds,"lanczos", "one", hs[i], 0);
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
    end
    println("ms: ", ms);

    end
    #TODO: plot it
    plot(hs, ms)
    savefig("./szplotnew.png")
end

function calculatefidelity()

    println("Starting fid!!");
    @time begin
    N=4;
    J=1
    deltah=0.001;
    hs=generateHListUniform(J, 20);
    fids=Any[];
    bonds = bondListFrustrated(N)
    println("hs: ", hs);
    for i=1:length(hs)
        println("starting h: ", hs[i])
        temp = calculateEigensystemTransverse(N, J, hs[i], bonds,"lanczos", "one", hs[i], 0);
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


function plottest()
    x=[0.1, 0.10069555500567187, 0.10139594797900289, 0.10210121257071929, 0.1028113826656066, 0.10352649238413769, 0.10424657608411206, 0.10497166836230663, 0.10570180405613792, 0.10643701824533586, 0.10717734625362917, 0.10792282365044256, 0.10867348625260563, 0.10942937012607376, 0.11019051158766086, 0.11095694720678427, 0.11172871380722174, 0.11250584846888068, 0.11328838852957958, 0.11407637158684206]
 y=[-0.8134175639031954, -0.8200063867643509, -0.8266616508941389, -0.8333842846052618, -0.8401752333849276, -0.8470354602567177, -0.8539659461491816, -0.8609676902725738, -0.8680417105021505, -0.8751890437692833, -0.8824107464601441, -0.8897078948215356, -0.8970815853746135, -0.9045329353360714, -0.9120630830466337, -0.9196731884078895, -0.9273644333256733, -0.9351380221620913, -0.9429951821939834, -0.9509371640799641]
    display(plot(x, y, seriestype=:scatter));
    savefig("./szplot.png");
end


function fidTest()
    println("FID TEST STARTING");
    println("FID TEST STARTING");
    N=4;
    J=1
    deltah=0.01;
    h=0.1;
    fids=Any[];
    bonds = bondListFrustrated(N)

    @time begin
        temp = calculateEigensystemTransverse(N, J, h, bonds,"lanczos", "one", h, 0);
        eigenvalues = temp[1]
        eigenvectors = temp[2]
        println("length of temp[3] ", length(temp[3]))
        println("length of temp ", length(temp))

        eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
        println("length of eigensystem ", length(eigensystem))
        println("eigensystem[3] ", eigensystem[3]);
        fid=calculateFidelity(eigensystem[2], temp[3][eigensystem[3]], h, N, deltah, bonds, J)
    end
    println("fid: ", fid);
end


function innerproducttest()
    println("Starting inner");
    eigenvector=[1/sqrt(2), 1/sqrt(2), 0, 0];
    eigenvector2=[1/sqrt(2), 1/sqrt(2), 0, 0];

    states=[1, 2, 3, 4];
    states2=[1, 3, 5, 6];
    println(innerProduct(eigenvector, states, eigenvector2, states2));
end
