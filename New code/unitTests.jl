include("entanglement helper.jl")


function test()
    N=2;
    bonds = bondListTwo();
    println(bonds);
    refStatesData=referenceStatesXY(N);

    println("TE SIZEH", length(refStatesData[1]));
    #list of all viable reference states
    refStates=refStatesData[1];
    #THIS maps a number to its state object, which contains info about ref state and shit
    refStatesMap=refStatesData[2];
    pt=momentum(0,0);
    viableSt=getViableStates2d(pt.px, pt.py, N, refStates);
    println("NEW");
    temp = calculateEigensystemHeisenbergMomentum2d(N, 1, bonds, "full", "all")
    println(temp[1]);
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

        hs=generateHListLog(J, 20);
        listA=plaquetteIndicies(generateCheckerboardPlaquettes(N)[1], N);
        bonds = bondListFrustrated(N)
        for i=1:length(hs)
            temp = calculateEigensystemTransverse(N, J, hs[i], bonds,"lanczos", "one", hs[i], 0);
            eigenvalues = temp[1]
            eigenvectors = temp[2]
            entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
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


function calculatesz()
    println("Starting sz!!");
    @time begin
    N=4;
    J=1
    hs=generateHListLog(J, 20);
    ms=Any[];
    bonds = bondListFrustrated(N)
    for i=1:length(hs)
        println("starting h: ", hs[i])
        temp = calculateEigensystemTransverse(N, J, hs[i], bonds,"lanczos", "one", hs[i], 0);
        eigenvalues = temp[1]
        eigenvectors = temp[2]
        sz=calculateSz(eigenvectors[1], temp[3][1], N);
        #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
        push!(ms, sz);
    end
    println("ms: ", ms);

    end
    #TODO: plot it
    plot(hs, sz)
    savefig("./szplot.png")
end
