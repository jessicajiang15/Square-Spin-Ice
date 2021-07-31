
include("unitTests.jl")
using Plots


function main()
    @time begin
        Random.seed!(123);
    println("hello")
    latticeType = "none"
    #heisenberg or transverse
    hamiltonianType = "transverse"
    #symmetry or momentum2d or reflection
    method = "none"
    #lanczos (Krylovit), full (LinearALgebra), sparse (Arpack)
    eigmethod="lanczos"
    #test=16, all=all, one=1, ignore if eigmethod=full
    num="one"
    file=false;

    bonds::Array{bond} = bond[]
    local eigenvalues
    local eigenvectors
    local temp
    J = 1;
    N = 4;
    h = 1;
    numSites=1;
    width = 0;



    if (latticeType == "frustrated")
        bonds = bondListFrustrated(N)
    elseif (latticeType == "four")
        bonds = bondListFourNeighbors(N)
    elseif (latticeType == "eight")
        bonds = bondListEightNeighbors(N)
    elseif( latticeType == "two")
        bonds=bondListTwo();
    end

    if (hamiltonianType == "heisenberg")
        if (method == "symmetry")
            temp = calculateEigensystemHeisenberg(N, J, bonds, eigmethod, num)
        elseif(method=="momentum")
            temp = calculateEigensystemHeisenbergMomentum(N, J, bonds, eigmethod, num)
        elseif(method=="momentum2d")
            temp = calculateEigensystemHeisenbergMomentum2d(N, J, bonds, eigmethod, num)
        elseif(method=="reflection")
            temp=calculateEigensystemHeisenbergReflection(N, J, bonds, eigmethod, num);
        else
            println("error");
        end
    elseif (hamiltonianType == "transverse")
        if(method=="none")
            temp=calculateEigensystemTransverseNoSymmetry(numSites, J, h, bonds,eigmethod, num, h, width, "H2");
        else
            temp = calculateEigensystemTransverse(numSites, J, h, bonds,eigmethod, num, h, width);
        end
    end

    eigenvalues = temp[1]
    eigenvectors = temp[2]
    listA=[0, 1, 4, 5];

    println("eigenvalues: ", (eigenvalues));
    println((eigenvectors));

    if(file)
        io=open("eigenvalues"*latticeType*hamiltonianType*method*eigmethod*".txt", "w") do io
            for i=1:length(eigenvalues)
                write(io, "$(eigenvalues[i]) \n")
            end
        end

        io=open("eigenvectors"*latticeType*hamiltonianType*method*eigmethod*".txt", "w") do io
            for i=1:length(eigenvectors)
                write(io, "``$(eigenvectors[i]) \n")
            end
        end
    end

    println("bye")
    println("TOTAL TIME");

end
end

println()

#main();
#fileiotest()
#thetest();
#println("starting");
#lol()
#entanglementtest();
#testbond();
#test();
#calculategsentanglement()
#println(singleOutEvenOddSpins(true,2^(4*4), 4))[1];
#thesztest()
#thesztest2()
#plaquetteTests()
#calculateStaggeredFlippabilityTest()

#fidTest2()
#plottest()
#calculateSpiTestAbs()
#calculateSpiTest()
#test();
#fidTest()
#fidTest2()
#plottest()
#testMeasures();

#testPlot();
#calculatefidelity()
#testExtractDigits()

#calculategsentanglement()


#graphTest();
#calculatesz()

#plottest()

#fidTest()

#innerproducttest()
#graphTest()
#thesztest2()
#fileiotest()

#filetest()
testlol();
#thesztestinfinitelattice()
