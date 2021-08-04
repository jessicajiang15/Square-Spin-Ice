include("entanglement helper.jl")
include("nlchelpers.jl")


function test()
    plaquettes=generateCheckerboardNoCrossPlaquettes(4);
    println(length(plaquettes));
    println(plaquettes);

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
    N=4
    eigenvector=Float64[1/sqrt(2), 0,0,0,0,0,0,0,0,0,0,0,0,0,0,1/sqrt(2)];
    println(length(eigenvector));
    listA=Int[1, 2];
    states::Array{Int}=Int[];
    for i=0:2^(N)-1
        push!(states, i);
    end
    println("calculating entropy");
    #eigenvector, states, listA, N
    println(getEntanglementEntropy(eigenvector, states, listA, N));
end

function lol()
        println("starting");
    println("lmao");
    N=7;
    listA=Int[5, 6, 7, 2];
    state=9
    #=
    N=2;
    listA=Int[0, 1];
    for i=0:15
        println(getBNum(i, listA, 2))
    end
=#

        getA=getANum(state, listA);
        getB=getBNum(state, listA, 4);
        println("a", getA);
        println("b", getB);
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

        hs=generateHListUniform(J, 50);
        listA=plaquetteIndicies(generateCheckerboardPlaquettes(N)[1], N);
        println(listA);
        bonds = bondListFrustrated(N)
        for i=1:length(hs)
            temp = calculateEigensystemTransverse(N*N, J, hs[i], bonds,"lanczos", "one", hs[i], 0);
            eigenvalues = temp[1]
            eigenvectors = temp[2]
            eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
            entropy=getEntanglementEntropy(eigensystem[2], temp[3][eigensystem[3]], listA, N*N);
            push!(entropies, entropy);
        end
    #println(entropy);
    end
    println(entropies)
    plot(hs, entropies)
    #savefig("Users/Jessica/git/Square Spin Ice/szplot.png")
    savefig("./entanglementplotnew.png")

end


function graphTest()
    x=[0.1, 0.14500000000000002, 0.19000000000000003, 0.23500000000000004, 0.28, 0.325, 0.37, 0.415, 0.45999999999999996, 0.505, 0.55, 0.5950000000000001, 0.6400000000000001, 0.6850000000000002, 0.7300000000000002, 0.7750000000000002, 0.8200000000000003, 0.8650000000000003, 0.9100000000000004, 0.9550000000000004]
y = [987.6543209876542, 987.6543209876542, 12.728366825097414, 987.6543209876542, 987.6543209876542, 987.6543209876542, 5.877008984993465, 987.6543209876542, 2.5292069246633604, 987.6543209876542, 987.6543209876542, 0.9811808433409414, 0.7661818838389369, 0.6130510365310856, 0.5004577612292375, 0.4153976405436494, 0.34963484689566593, 987.6543209876542, 987.6543209876542, 0.2222886060269335]
plot(x, y)
display(plot(x, y))
#savefig("./plot.png")
end




function plottest()
    x=[0.1, 0.12000000000000001, 0.14, 0.16, 0.18, 0.19999999999999998, 0.21999999999999997, 0.23999999999999996, 0.25999999999999995, 0.27999999999999997, 0.3, 0.32, 0.34, 0.36000000000000004, 0.38000000000000006, 0.4000000000000001, 0.4200000000000001, 0.4400000000000001, 0.46000000000000013, 0.48000000000000015, 0.5000000000000001]


 y=[0.025435207891312483, 0.025484529594182126, 0.03693435913257226, 0.06128172661284137, 0.09468594876745073, 0.13182056706838574, 0.16846369415488527, 0.20218476528790308, 0.23200479713697475, 0.25782671211855523, 0.27998297083289536, 0.29895999632637915, 0.3152568114936234, 0.3293241610033619, 0.3415455506998788, 0.3522375955642083, 0.36165788615940253, 0.3700147716144838, 0.3774766718516003, 0.38418008097641, 0.39023611757754706]
    display(plot(x, y));
    savefig("./szplotinflattice, order 20.png");
end



function innerproducttest()
    println("Starting inner");
    eigenvector=[1/sqrt(2), 1/sqrt(2), 0, 0];
    eigenvector2=[1/sqrt(2), 1/sqrt(2), 0, 0];

    states=[1, 2, 3, 4];
    states2=[1, 3, 5, 6];
    println(innerProduct(eigenvector, states, eigenvector2, states2));
end



function thesztest2()
    println("Starting sz!!");
    @time begin
    N=16;
    J=1
    h=0.25
    ms=Any[];
    bonds = bondListFrustrated(4)

        temp = calculateEigensystemTransverse(N, J, h, bonds,"lanczos", "one", h, 0);
        eigenvalues = temp[1]
        eigenvectors = temp[2]
        eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);

        println("length of temp[3] ", length(temp[3]))
        println("length of temp ", length(temp))
        println("length of eigensystem ", length(eigensystem))
        println("eigensystem[3] ", eigensystem[3]);
        sz=calculateSz(eigensystem[2], temp[3][eigensystem[3]], N);
        szm=szMatrix(eigensystem[2], temp[3][eigensystem[3]], N);
        #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);

    println("ms: ", sz);
    println("ms matrix: ", szm);
end
end

function fidTest2()
    println("starting fid");
    println("starting fid");


    num=50;
    N=4;
    J=1;
    hmin=0.1*J;
    hmax=J;
    bonds = bondListFrustrated(N)
    hs=generateHListUniform(hmin, hmax, num);
    println("hs", hs);

    fids=calculateFidelity(hmin, hmax, num, N, J, bonds);
    println("fids: ", fids);
    plot(hs, fids);
    savefig("./fidelityplot2.png");
end



function sPi()
    println("Starting sz!!");
    @time begin
    N=4;
    J=1
    h=0.35
    ms=Any[];
    bonds = bondListFrustrated(N)
    spis=Any[];
temp=calculateEigensystemTransverseNoSymmetry(N, J, h, bonds,"lanczos", "one", h, width, "H1");
        #temp = calculateEigensystemTransverse(N, J, h, bonds,"lanczos", "one", h, 0);
        eigenvalues = temp[1]
        eigenvectors = temp[2]
        eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);

        println("length of temp[3] ", length(temp[3]))
        println("length of temp ", length(temp))
        println("length of eigensystem ", length(eigensystem))
        println("eigensystem[3] ", eigensystem[3]);
        spi=calculateSPiSz(eigensystem[2], temp[3][eigensystem[3]], N, temp[4][eigensystem[3]]);
        push!(spis, spi);
        #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);

    println("ms: ", sz);
    println("ms matrix: ", szm);
end
end





function calculateSpiTestAbs()
    println("Starting sz!!");
    println("Starting sz!!");

    @time begin

    N=4;
    J=1
    hs=generateHListUniform(0.1, 1, 50);
    spis=Any[];
    bonds = bondListFrustrated(N)
    println("hs: ", hs);
    for i=1:length(hs)
        println("starting h: ", hs[i])
        temp =calculateEigensystemTransverseNoSymmetry(N, J, hs[i], bonds,"lanczos", "one", hs[i], 0, "H1");
        eigenvalues = temp[1]
        eigenvectors = temp[2]
        #eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);

        #println("length of temp[3] ", length(temp[3]))
        #println("length of temp ", length(temp))
        #println("length of eigensystem ", length(eigensystem))

        #println("eigensystem[3] ", eigensystem[3]);
        println("starting time");
        @time begin
        spi=calculateSPiSzNewAbs(eigenvectors[1], temp[3][1], N);
    end
        println("h: ", hs[i], " , spi: ", spi);

        push!(spis, spi);

    end
    println("spis: ", spis);
    end
    #TODO: plot it
    plot(hs, spis)
    savefig("./spiAbsplot.png")
end


function plaquetteTests()
    N=4
    squareIndicies=generateListsofPlaquetteIndiciesFlip(N);
    println("Square indicies: ", squareIndicies);
    state=9508;
    for i=1:length(squareIndicies)
        println("plaquette: ", squareIndicies[i], " , number: ", getPlaquetteNumber(squareIndicies[i][1], N))
        println("is neel: ", isNeel(state, squareIndicies[i], N));
    end
end

function fileiotest()
r=readdlm("graphs.txt"; skipblanks=false);
println(r[10, 1]);
#println(r);
#println((r));
#println((r));
#=
count=1;
num=r[1, count];
while(num!="")
    if(count==1)
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
    num=r[1, count];
    println(num);
    count+=1;
end
        #close(file);
        =#
end


function filetest()
    d=readdlm("graphs.txt"; skipblanks=false);
    r=readFromGraphFile();
    println(r);
    #println(r);
end

function testlol()
    @time begin
    order=56;
    graphs=readFromGraphFile();
    sz=calculateInfiniteLatticeSz(order, 1, 10, graphs, 0)
    println("sz ", sz);
end
end


function huh()
    N=4;
    listA=plaquetteIndicies(generateCheckerboardPlaquettes(N)[1], N);
    println(listA);
end



function entanglementtest2()
    println("Starting!!");
    @time begin
    N=4;
    J=1
    J2=1;
    h=0.1

        listA=plaquetteIndicies(generateCheckerboardNoCrossPlaquettes(N)[1], N);
        println(listA);
        bonds = bondListFrustrated(N)
            temp = calculateEigensystemTransverse(N*N, J, J2, h, bonds,"lanczos", "one", h, 0);
            eigenvalues = temp[1]
            eigenvectors = temp[2]
            eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
            println(eigenvalues);
            entropy=getEntanglementEntropy(eigensystem[2], temp[3][eigensystem[3]], listA, N*N);
    #println(entropy);
    println(entropy)
    #savefig("Users/Jessica/git/Square Spin Ice/szplot.png")
end

end


function spitest()
    println("Starting sz!!");
    println("Starting sz!!");

    @time begin
    N=4;
    J=1
    J2=1;
    spis=Any[];
    h=0.1
    bonds = bondListFrustrated(N)
        temp =calculateEigensystemTransverseNoSymmetry(N*N, J, J2, h, bonds,"lanczos", "one", h, 0, "H1");
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
        println("spi", spi);
    end
    end
end
