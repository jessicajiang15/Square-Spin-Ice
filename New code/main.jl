
include("unitTests.jl")
using Plots
#ENV["JULIA_NUM_THREADS"] = 12

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

ctc_calculation()
#plot_transport_norm_distribution()
#plot_distribution_localization_length()
#plot_transport_norm_distribution()
#plot_transport_norm_distribution()
#plot_transport_norm_vs_L_quasiperiodic_sampling()
#plot_O_distances_distribution()
#plot_transport_norm_vs_L_quasiperiodic()
#plot_transport_norm_vs_L_quasiperiodic_sampling()
#plot_Ws_distances_distribution()
#plot_Ws_distances_distribution()
#plot_transport_norm_vs_L_quasiperiodic_sampling()
#time_test_interacting()
#plot_transport_norm_vs_i0()
#plot_energy_distributions()
#plot_v_distribution()
#plot_v_distribution()
#plot_orbital_correction_norm_avg()
#plot_transport_norm_vs_L_quasiperiodic()
#plot_transport_norm_vs_L_quasiperiodic()
#plot_from_data_mean()
#plot_distribution_localization_length()
#plot_distribution_localization_length()
#plot_energy_distribution()
#print(to_sum_k_rec(2, 10))
#anderson_diagonalization()
#anderson_diagonalization_new()
#test_energy_matrix()
#plot_from_data_mean()
#w_svd_test()
#plot_w_norm_from_data()
#print_Ws_max()
#plot_distribution_localization_length()
#plot_from_data()
#plot_orbital_correction_norm_avg()
#plot_from_data_2()
#plot_O_distances_distribution()
#quasiperiodic_diagonalization_transition()
#plot_from_data()
#plot_transport_norm_vs_L_disorder()
#plot_transport_norm_vs_L_quasiperiodic()
#plot_from_data_mean()
#plot_transport_norm_vs_L_quasiperiodic()
#plot_transport_norm_average_distribution()
#plot_transport_norm_vs_V_quasiperiodic()
#plot_from_data()
#plot_transport_norm_vs_L_disorder()22
#plot_from_data()
#plot_transport_norm_vs_L_quasiperiodic()
#plot_transport_norm_vs_L_quasiperiodic()
#plot_transport_norm_vs_V_disorder()
#plot_transport_norm_average_distribution()
#plot_transport_norm_distribution()
#plot_Ws_distances_distribution()
#plot_g_norms_as_function_distance()
#plot_g_norms_as_function_distance()
#plot_quasiperiodic_norm_distribution()
#quasiperiodic_diagonalization()
#plot_Ws_norm_distribution()
#plot_from_data_avglog()
#plot_w_norm_from_data()
#plot_Ws_norm_distribution()
#plot_Ws_norm_distribution()
#plot_Ws_distribution()
#time_Ws()
#plot_from_data_avglog()
#plot_w_norm_from_data()
#plot_distribution_deroeck()
#plot_Ws_distribution()
#anderson_diagonalization()
#plot_distribution_deroeck()
#plot_Fs_norm_distribution()
#plot_from_data()
#plot_quasiperiodic_norm_distribution()
#plot_Ws_distribution_quasiperiodic()
#plot_Ws_norm_distribution()