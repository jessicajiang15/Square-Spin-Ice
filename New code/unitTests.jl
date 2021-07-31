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

        hs=generateHListUniform(J, 50);
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
    x=[0.1, 0.10069555500567187, 0.10139594797900289, 0.10210121257071929, 0.1028113826656066, 0.10352649238413769, 0.10424657608411206, 0.10497166836230663, 0.10570180405613792, 0.10643701824533586, 0.10717734625362917, 0.10792282365044256, 0.10867348625260563, 0.10942937012607376, 0.11019051158766086, 0.11095694720678427, 0.11172871380722174, 0.11250584846888068, 0.11328838852957958, 0.11407637158684206]


 y=[6.739980572365193, 7.908802091581995, 9.3482730178918, 11.039936762469443, 12.89465316346356, 14.72099147833389, 16.22635839736444, 17.08560667795857, 17.071859626950054, 16.172737502414584, 14.598062125143258, 12.671841414404305, 10.696958715571288, 8.876861035898184, 7.307936113547052, 6.009702345012287, 4.960041302269914, 4.1206974452076714, 3.451528254105747, 2.9168214844671394, 2.4871901079529612, 2.1394126568707184, 1.8555177293119485, 1.6217265647538384, 1.4274997861270966, 1.2647595864617651, 1.1272854653317008, 1.0102565317336325, 0.9099091333490904, 0.823282004444078, 0.7480266077108942, 0.6822655574725798, 0.6244864624867144]
 y1=[0.5734618301059907, 0.5281882387390504, 0.48783980805195054, 0.4517323616624759, 0.41929559306832886, 0.3900513352347353, 0.363596435143867, 0.33958919888113465, 0.3177385651302, 0.297795439837723, 0.2795456621438757, 0.26280433244342927, 0.2474111556306625, 0.23322661654963156, 0.2201288525555328, 0.20801103193410936, 0.19677917871569392];
 append!(y, y1);
    display(plot(x, y, seriestype=:scatter));
    savefig("./szplot.png");
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
    order=56;
    graphs=readFromGraphFile();
    sz=calculateInfiniteLatticeSz(order, 1, 10, graphs, 0)
    println("sz ", sz);
end
