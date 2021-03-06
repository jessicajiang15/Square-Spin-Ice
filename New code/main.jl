
include("unitTests.jl")
using Plots


function main()
    @time begin
        Random.seed!(123)
        println("hello")
        latticeType = "frustrated"
        #heisenberg or transverse
        hamiltonianType = "transverse"
        #symmetry or momentum2d or reflection
        method = "symmetry"
        #lanczos (Krylovit), full (LinearALgebra), sparse (Arpack)
        eigmethod = "lanczos"
        #test=16, all=all, one=1, ignore if eigmethod=full
        num = "one"
        file = false
        bonds::Array{bond} = bond[]
        local eigenvalues
        local eigenvectors
        local temp
        J = 1
        J2 = 1
        N = 4
        h = 1
        numSites = 16
        width = 0



        if (latticeType == "frustrated")
            bonds = bondListFrustrated(N)
        elseif (latticeType == "four")
            bonds = bondListFourNeighbors(N)
        elseif (latticeType == "eight")
            bonds = bondListEightNeighbors(N)
        elseif (latticeType == "two")
            bonds = bondListTwo()
        end
        println("BONDS LENGTH??", length(bonds))

        if (hamiltonianType == "heisenberg")
            if (method == "symmetry")
                temp =
                    calculateEigensystemHeisenberg(N, J, bonds, eigmethod, num)
            elseif (method == "momentum")
                temp = calculateEigensystemHeisenbergMomentum(
                    N,
                    J,
                    bonds,
                    eigmethod,
                    num,
                )
            elseif (method == "momentum2d")
                temp = calculateEigensystemHeisenbergMomentum2d(
                    N,
                    J,
                    bonds,
                    eigmethod,
                    num,
                )
            elseif (method == "reflection")
                temp = calculateEigensystemHeisenbergReflection(
                    N,
                    J,
                    bonds,
                    eigmethod,
                    num,
                )
            else
                println("error")
            end
        elseif (hamiltonianType == "transverse")
            if (method == "none")
                temp = calculateEigensystemTransverseNoSymmetry(
                    numSites,
                    J,
                    J2,
                    h,
                    bonds,
                    eigmethod,
                    num,
                    h,
                    width,
                    "H2",
                )
            else
                temp = calculateEigensystemTransverse(
                    numSites,
                    J,
                    J2,
                    h,
                    bonds,
                    eigmethod,
                    num,
                    h,
                    width,
                )
            end
        end

        eigenvalues = temp[1]
        eigenvectors = temp[2]
        listA = [0, 1, 4, 5]

        println("eigenvalues: ", (eigenvalues))
        #println(eigenvalues[1]==eigenvalues[2]);
        #println((eigenvectors));

        if (file)
            io = open(
                "eigenvalues" *
                latticeType *
                hamiltonianType *
                method *
                eigmethod *
                ".txt",
                "w",
            ) do io
                for i = 1:length(eigenvalues)
                    write(io, "$(eigenvalues[i]) \n")
                end
            end

            io = open(
                "eigenvectors" *
                latticeType *
                hamiltonianType *
                method *
                eigmethod *
                ".txt",
                "w",
            ) do io
                for i = 1:length(eigenvectors)
                    write(io, "$(eigenvectors[i]) \n")
                end
            end
        end

        println("bye")
        println("TOTAL TIME")

    end
end

println()
#ultimatetestttt()
#findHcFromIterations()
#trackEDMeanFieldSzEachIteration()
meanFieldSzInfiniteLatticeManyJ2AndOrders()
#lel();
#metropolisTest()
#bondsmaptest()
#meanFieldTest()
#szMeanFieldInfiniteLattice()
#testtest()
#trackEDMeanFieldSzEachIteration()
#calculateEDMeanFieldSz()
#testtest()
#plaquetteTest()
#bondsTest()
#meanFieldSzInfiniteLatticeManyJ2AndOrders()
#fidelityInfiniteLattice()
#susinfinitelatticemanyJ2andorders()
#fidelityInfiniteLatticeMultipleJ2()
#fidelityInfiniteLattice()
#scalez2()fidTest2()
#szinfinitelatticemanyJ2andorders()
#thesztestj1j2ED()
#scalesz()
#szinfinitelatticenew()
#susceptibilityJ2J1ED()
#entanglementEDj1j2vsh()
#szinfinitelatticemanyJ2andorders()
#entanglementinfinitelatticemultipleJ2andorders()
#spiinfinitelatticemultipleJsandorders()
#entanglementinfinitelatticenewj1j2()
#entanglementinfinitelatticenew()
#theentanglementtestinfinitelattice()
#testNumbering()
#main();
#thesztestinfinitelattice()
#spitest()
#spij1j2()
#weightspigraphs()
#spiinfinitelatticeone()
#spiinfinitelattice()
#spiinfinitelatticeold()
#sPi();
#flippabilityj1j2()
#testlo()
#spiinfinitelattice()
#plottest()
#entanglementasj1j2()
#spij1j2()
#entanglementtest()
#
#entanglementinfinitelatticenew()
#theentanglementtestinfinitelattice()
#entanglementtest2()
#entanglemenvttest();
#testbond();
#spiinfinitelattice()
#spiGraphj20()
#spiGraphj201()
#spiGraphj2028()
#spiGraphj2046()
#spiGraphj2082()
#spiGraphj21()
#spiinfinitelattice()
#test();
#spiNLCj1j2()
#calculateSpij1j2ED()
#calculategsentanglement()
#println(singleOutEvenOddSpins(true,2^(4*4), 4))[1];
#thesztest()
#thesztest2()
#plaquetteTests()
#szinfinitelatticewithED()
#calculateStaggeredFlippabilityTest()
#entanglementinfinitelatticenew()
#fidTest2()
#alculateSpij1j2ED()
#plottest()
#calculateSpiTestAbs()
#calculateSpiTest()
#test();
#fidTest()
#fidTest2()
#szNLCj1j2()
#plottest()
#testMeasures();
#sPi()
#testMeasures()
#testPlot();
#calculatefidelity()
#calculatefidelityJ1J2ED()
#fidTest2()
#testExtractDigits()
#lol()
#szinfinitelatticenew()
#calculateSpiTest()
#calculategsentanglement()
#entanglementinfinitelatticenew()
#entanglementinfinitelatticenew()
#graphTest();
#thesztest()
#calculatesz()

#plottest()

#fidTest()

#innerproducttest()
#graphTest()
#thesztest2()
#fileiotest()
#lol();
#filetest()
#weightsszgraphs()
#testlol();
#thesztestinfinitelattice()
#weightsszgraphs2()
#weightsszgraphs3()
#spiGraphj2064()
#calculateBaseWeightMeanFieldSz(1, 1, 1, 0.1, 10)
#spiinfinitelattice()
#spiinfinitelatticemultipleJsandorders()
#weightsszgraphs4()
#moreTest()
#spiinfinitelatticemultipleJsandorders()
#spiGraphj201()
#spiGraphj201()
#spiGraphj2082()
#spiGraphj20ratiotest()
#spiGraphj2064ratiotest()
#spiGraphj201ratiotest()
#spiGraphj2028ratiotest()
#spiGraphj2082ratiotest()
#moreTest()
#spiGraphj2064ratiotest()
#spiGraphj2046ratiotest()
#spiinfinitelattice()
#spiGraphj201()
#spiGraphj2082()
#hamiltonianTest()
