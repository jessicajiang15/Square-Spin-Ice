include("hamiltonianDiagonalization.jl")
include("calculations and plotting.jl")

function plot()

end

function main()
    @time begin
    println("hello")
    latticeType = "four"
    hamiltonianType = "heisenberg"
    method = "symmetry"

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
        elseif(method=="momentum")
            temp = calculateEigensystemHeisenbergMomentum(N, J, bonds)
        else
            temp = calculateEigensystemHeisenbergMomentum2d(N, J, bonds)
        end
    elseif (hamiltonianType == "transverse")
        temp = calculateEigensystemTransverse(N, J, h, bonds)
    end
    eigenvalues = temp[1]
    eigenvectors = temp[2]
    println((eigenvalues))
    #println((eigenvectors));
    println("bye")
    println("TOTAL TIME");
end
end

function test()
    println("NEW");
    N = 4
    temp = 6574;
    #println(rotateAllBits(1, temp, N));
    #referenceStates(N);
end
main();
#test();
