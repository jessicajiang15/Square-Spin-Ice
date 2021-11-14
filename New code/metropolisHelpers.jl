using Random
include("necessaryStructures.jl")
include("necessaryBitOperations.jl")


function randomConfiguration(N)
    return rand(0:2^N-1);
end

function calculateIsingEnergy(state, bonds, J)
    energy=0;
    for i=1:length(bonds)
        site1=bonds[i].site1;
        site2=bonds[i].site2;
        digit1=getTi(site1.num, state);
        a1= digit1==1 ? 1/2 : -1/2;
        digit2=getTi(site2.num, state);
        a2= digit2==1 ? 1/2 : -1/2;
        spinCoupling=J*a1*a2;
        energy+=spinCoupling;
    end
    return energy;
end

function metropolisAlgorithm(thermalizationSteps, iterations, bonds, N, J, T)
    state=randomConfiguration(N);
    energy=calculateIsingEnergy(state, bonds, J)
    energies=Float64[];
    for i=1:iterations
        if(i>=thermalizationSteps)
            push!(energies, energy);
        end
        for i=0:N-1
            b::Int=flipBit(i, state);
            temp=calculateIsingEnergy(b, bonds, J);
            if(temp<energy)
                energy=temp;
                state=b;
            else
                probability=exp((1/T)*energy)*exp((1/T)*temp);
                if(rand()<probability)
                    energy=temp;
                    state=b;
                end
            end
        end
    end
    return energies;
end
