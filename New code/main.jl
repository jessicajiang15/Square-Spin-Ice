include("hamiltonianDiagonalization.jl")
include("calculations and plotting.jl")

function plot()

end

function main()
    @time begin
        Random.seed!(123);
    println("hello")
    latticeType = "four"
    hamiltonianType = "transverse"
    method = "symmetry"
    #lanczos (Krylovit), full (LinearALgebra), sparse (Arpack)
    eigmethod="full"
    #test=16, all=all, one=1, ignore if eigmethod=full
    num="all"
    file=true;

    bonds::Array{bond} = bond[]
    local eigenvalues
    local eigenvectors
    local temp
    J = 1;
    N = 4;
    h = 1;
    width = 0;


    if (latticeType == "frustrated")
        bonds = bondListFrustrated(N)
    elseif (latticeType == "four")
        bonds = bondListFourNeighbors(N)
    elseif (latticeType == "eight")
        bonds = bondListEightNeighbors(N)
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
        temp = calculateEigensystemTransverse(N, J, h, bonds,eigmethod, num, h, width);
    end

    eigenvalues = temp[1]
    #eigenvectors = temp[2]

    println((eigenvalues));
    #println((eigenvectors));
    if(file)
        io=open("eigenvalues"*latticeType*hamiltonianType*method*".txt", "w") do io
            for i=1:length(eigenvalues)
                write(io, "$(eigenvalues[i]) \n")
            end
        end

        io=open("eigenvectors"*latticeType*hamiltonianType*method*".txt", "w") do io
            for i=1:length(eigenvectors)
                write(io, "$(eigenvectors[i]) \n")
            end
        end
    end

    println("bye")
    println("TOTAL TIME");

end
end

function test()
    println("NEW");
    d=Normal(1, 0);
    println(rand(d, 10));
end

main();

#test();
