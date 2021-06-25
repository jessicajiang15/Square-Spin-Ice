include("hamiltonianDiagonalization.jl")
include("calculations and plotting.jl")

function plot()

end

function main()
    println("hello")
    latticeType = "frustrated"
    hamiltonianType = "heisenberg"
    method = "momentum"

    bonds::Array{bond} = bond[]
    local eigenvalues
    local eigenvectors
    local temp
    J = 1;
    N = 4;
    h = 1;

    if (latticeType == "frustrated")
        bonds = bondListFrustrated(N)
    elseif (latticeType == "four")
        bonds = bondListFourNeighbors(N)
    elseif (latticeType = "eight")
        bonds = bondListEightNeighbors(N)
    end

    if (hamiltonianType == "heisenberg")
        if (method == "symmetry")
            temp = calculateEigensystemHeisenberg(N, J, bonds)
        else
            temp = calculateEigensystemHeisenbergMomentum(N, J, bonds)
        end
    elseif (hamiltonianType == "transverse")
        temp = calculateEigensystemTransverse(N, J, h, bonds)
    end
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    println((eigenvalues))
    #println((eigenvectors));

    println("bye")


end

function test()
    println("NEW");
    N = 4
    temp = 6574;
    println(rotateAllBits(1, temp, N));
end
main();
#test();
