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
    J2=0
    entropies=Any[];
    hmin=0.1;
    hmax=2;

        hs=generateHListUniform(hmin, hmax, 100);
        listA=plaquetteIndicies(generateCheckerboardPlaquettes(N)[1], N);
        println(listA);
        bonds = bondListFrustrated(N)
        for i=1:length(hs)
            temp = calculateEigensystemTransverse(N*N, J, J2, hs[i], bonds,"lanczos", "one", hs[i], 0);
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
    savefig("./entanglementplot J2: "*string(J2)*", hmin: "*string(hmin)*", hmax: "*string(hmax)*".png");

end


function graphTest()
    x=[0.1, 0.14500000000000002, 0.19000000000000003, 0.23500000000000004, 0.28, 0.325, 0.37, 0.415, 0.45999999999999996, 0.505, 0.55, 0.5950000000000001, 0.6400000000000001, 0.6850000000000002, 0.7300000000000002, 0.7750000000000002, 0.8200000000000003, 0.8650000000000003, 0.9100000000000004, 0.9550000000000004]
y = [987.6543209876542, 987.6543209876542, 12.728366825097414, 987.6543209876542, 987.6543209876542, 987.6543209876542, 5.877008984993465, 987.6543209876542, 2.5292069246633604, 987.6543209876542, 987.6543209876542, 0.9811808433409414, 0.7661818838389369, 0.6130510365310856, 0.5004577612292375, 0.4153976405436494, 0.34963484689566593, 987.6543209876542, 987.6543209876542, 0.2222886060269335]
plot(x, y)
display(plot(x, y))
#savefig("./plot.png")
end




function plottest()
    x=[0.1, 0.138, 0.17600000000000002, 0.21400000000000002, 0.252, 0.29, 0.32799999999999996, 0.36599999999999994, 0.4039999999999999, 0.4419999999999999, 0.47999999999999987, 0.5179999999999999, 0.5559999999999999, 0.594, 0.632, 0.67, 0.7080000000000001, 0.7460000000000001, 0.7840000000000001, 0.8220000000000002, 0.8600000000000002, 0.8980000000000002, 0.9360000000000003, 0.9740000000000003, 1, 1.0120000000000002, 1.0500000000000003, 1.0880000000000003, 1.1260000000000003, 1.1640000000000004, 1.2020000000000004, 1.2400000000000004, 1.2780000000000005]
    x1=[1.3160000000000005, 1.3540000000000005, 1.3920000000000006, 1.4300000000000006, 1.4680000000000006, 1.5060000000000007, 1.5440000000000007, 1.5820000000000007, 1.6200000000000008, 1.6580000000000008, 1.6960000000000008, 1.7340000000000009, 1.772000000000001, 1.810000000000001, 1.848000000000001, 1.886000000000001, 1.924000000000001, 1.962000000000001]

    append!(x, x1);

    println(length(x))

 y=[0.6931554973559289, 0.6931566513135997, 0.6931580030135706, 0.6931595917017824, 0.6931614655626608, 0.6931636840503151, 0.6931663209145511, 0.6931694681579544, 0.6931732412518335, 0.6931777860709394, 0.693183288199999, 0.6931899855520833, 0.6931981856727243, 0.6932082897741164, 0.6932208266069516, 0.693236501021826, 0.6932562650702393, 0.6932814250104411, 0.6933138088521826, 0.693356046260912, 0.6934120973550814, 0.6934885492062778, 0.6936001527141216, 0.693848810805625, 1.7896261747186408, 1.3881081056276219, 1.3863462082902736, 1.3862597575483824, 1.3862399529019687, 1.3862350147667382, 1.3862354160773416, 1.386238100437923, 1.3862417742713433]
 y1=[1.3862458065457506, 1.3862498686830083, 1.3862537872721254, 1.3862574741067692, 1.386260889496157, 1.3862640217425914, 1.3862668751786122, 1.386269463017255, 1.386271803024712, 1.3862739148966676, 1.3862758186831787, 1.3862775338666122, 1.3862790788484831, 1.386280470692036, 1.3862817250233803, 1.3862828560291103, 1.386283876510647, 1.386284797969876]

 append!(y, y1);
 println(length(y))
    display(plot(x, y));
    savefig("./j1j2entanglement");
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
    h=0.1
    ms=Any[];
    bonds = bondListFrustrated(4)

        temp = calculateEigensystemTransverse(N, J, 1, h, bonds,"lanczos", "one", h, 0);
        eigenvalues = temp[1]
        eigenvectors = temp[2]
        eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);

        println("length of temp[3] ", length(temp[3]))
        println("length of temp ", length(temp))
        println("length of eigensystem ", length(eigensystem))
        println("eigensystem[3] ", eigensystem[3]);
        sz=calculateSz(eigensystem[2], temp[3][eigensystem[3]], N);
        #szm=szMatrix(eigensystem[2], temp[3][eigensystem[3]], N);
        #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);

    println("ms: ", sz);
    #println("ms matrix: ", szm);
end
end

function fidTest2()
    println("starting fid");
    println("starting fid");


    num=50;
    N=4;
    J=1;
    J2=0;
    hmin=0.2;
    hmax=1;
    bonds = bondListFrustrated(N)
    hs=generateHListUniform(hmin, hmax, num);
    println("hs length: ", length(hs));

    println("hs", hs);

    fids=calculateFidelity(hmin, hmax, num, N, J, J2, bonds);
    println("fids: ", fids);
    println("fids length: ", length(fids));

    plot(hs, fids);
    savefig("./fidelityplot J2 is 0.png");
end



function sPi()
    println("Starting sz!!");
    @time begin
    N=4;
    J=1
    J2=1
    h=0.1
    bonds = bondListFrustrated(N)
    spis=Any[];
    temp =calculateEigensystemTransverseNoSymmetry(N*N, J, J2, h, bonds,"lanczos", "one", h, 0, "H1");
        #temp = calculateEigensystemTransverse(N, J, h, bonds,"lanczos", "one", h, 0);
        eigenvalues = temp[1]
        eigenvectors = temp[2]
        eigensystem=getLowestLyingStates(eigenvalues, eigenvectors);
        println("eigenvalues ", eigenvalues);
        spi=calculateSPiSzNew(eigenvectors[1], temp[3][1], N^2);
        #entropy=getEntanglementEntropy(eigenvectors[1], temp[3][1], listA, N);
        println("sipi ", spi);
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


function testNumbering()
    for i=1:16
        println("i, ", i, ", num: ", getPlaquetteNumber(i, 4));
    end
end

function testlo()
    graphs=readFromGraphFile();

    println(getLastGraphNumOrder(6, graphs));
end
