include("hamiltonianDiagonalization.jl")
include("calculations and plotting.jl")

function plot()

end

function main()
    @time begin
        Random.seed!(123);
    println("hello")
    latticeType = "four"
    #heisenberg or transverse
    hamiltonianType = "heisenberg"
    #symmetry or momentum2d or reflection
    method = "momentum2d"
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
        temp = calculateEigensystemTransverse(N, J, h, bonds,eigmethod, num, h, width);
    end

    eigenvalues = temp[1]
    eigenvectors = temp[2]

    println((eigenvalues));
    #println((eigenvectors));
    if(file)
        io=open("eigenvalues"*latticeType*hamiltonianType*method*eigmethod*".txt", "w") do io
            for i=1:length(eigenvalues)
                write(io, "$(eigenvalues[i]) \n")
            end
        end

        io=open("eigenvectors"*latticeType*hamiltonianType*method*eigmethod*".txt", "w") do io
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
    N=2;
    bonds = bondListTwo();
    println(bonds);
    refStatesData=referenceStatesXY(N);

    println("TE SIZEH", length(refStatesData[1]));
    #list of all viable reference states
    refStates=refStatesData[1];
    #THIS maps a number to its state object, which contains info about ref state and shit
    refStatesMap=refStatesData[2];
    pt=momentum(0,0);
    viableSt=getViableStates2d(pt.px, pt.py, N, refStates);
    println("NEW");
    temp = calculateEigensystemHeisenbergMomentum2d(N, 1, bonds, "full", "all")
    println(temp[1]);
end

main();

#test();
