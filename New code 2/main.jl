
include("unitTests.jl")
using Plots


function main()
    @time begin
        Random.seed!(123)
        println("hello")
        latticeType = "1d"
        #heisenberg or transverse
        hamiltonianType = "transverse"
        #symmetry or momentum2d or reflection
        method = "none"
        #lanczos (Krylovit), full (LinearALgebra), sparse (Arpack)
        eigmethod = "full"
        #test=16, all=all, one=1, ignore if eigmethod=full
        num = "all"
        file = false
        pbc::Bool=false
        bonds::Array{bond} = bond[]
        hs = Float64[]
        hs2 = Float64[]
        hss=generateHListUniformIncludeOne(0.1, 1, 10);
        print(hss)

        local eigenvalues
        local eigenvectors
        local temp
        J = 1
        J2 = 0
        N = 2
        h = 1
        numSites = 2
        width = 0

        for i=1:N
            push!(hs,h)
            push!(hs2,0)
        end


        if (latticeType == "frustrated")
            bonds = bondListFrustrated(N)
        elseif (latticeType == "four")
            bonds = bondListFourNeighbors(N)
        elseif (latticeType == "eight")
            bonds = bondListEightNeighbors(N)
        elseif (latticeType == "two")
            bonds = bondListTwo()
        elseif (latticeType=="1d")
            bonds=bonds1D(N, pbc)
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
                #(N, J, J2, hs::Vector{Float64},hs2::Vector{Float64},
                # bonds,eigmethod, num, h1orh2)
                temp = calculateEigensystemTransverseNoSymmetry(
                    numSites,
                    J,
                    J2,
                    hs,
                    hs2,
                    bonds,
                    eigmethod,
                    num,
                    "H1",
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

        println(temp)

        eigenvalues = temp[1]
        eigenvectors = temp[2]
        #println("eigenvalues: ", (eigenvalues))
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
#find_1d()
#find_1d_gs_per_site_43()
#find_1d_gs_per_site_diff_43()
#find_energies_121()
#fitttttt()
#find_diff_energy_44()
#calculate_correlations_121()
find_symmetry_energies_splitting_121()
#find_symmetry_energies_121()
#calculate_correlations_121()
#calculate_correlations_121()
#find_1d_gs_per_site_43()
#find_sz_121()
#calculateFidelity121()
