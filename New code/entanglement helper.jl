include("eigenvectorQuantityPlots.jl")

function getANum(state, listA)
    num=0;
    for i=1:length(listA)
        digit=getTi(listA[i]-1, state);
        num=num|(digit<<(i-1));
    end
    #println(num);
    return num;
end

#e.g. 4, 1 gives 00...0011110
function createMask(digits, pos)
    if(digits<=0)
        #println("wtf");
        return 0;
    end
    num=(2^(digits))-1;
    return num<<pos;
end

#create a number by taking out the bits at positions in listA
function getBNum(state, listA, N)
    #idea: have a mask for each and keep track of last digit that youve finished.
#shift to that digit
#println("state", state);
    num=0;
    lastNum=-1;
    digitsdone=0;
    for i=1:length(listA)
        temp=0;
        numDigits=listA[i]-lastNum-2;
        #println("numdigits",numDigits, "i, ", listA[i]);
        if(numDigits>0)
            mask=createMask(numDigits, lastNum+1);
            #println("mask",mask);
            #println("lastNum",lastNum);
            temp=(mask&state)>>(lastNum+1);
            #println("temp",temp);
            num|=temp<<digitsdone;
            #println("num",temp);
            digitsdone+=numDigits;
            #println("digits done ",digitsdone);
        end
        lastNum=listA[i]-1;
        #println("");
    end
    if(N*N-1!=listA[length(listA)-1])
        numDigits=N*N-lastNum-1;
        #println("numDigits",numDigits);
        if(numDigits==0)
            return num;
        end
        mask=createMask(numDigits, lastNum+1);
        temp=(mask&state)>>(lastNum+1);
        #println("mask&state",mask&state);
        #println("lastNum+1",lastNum+1);
        num|=temp<<digitsdone;
    end
    #println("final:", num);
    return num;
end


#create a number by taking out the bits at positions in listA
function getBNumStateNum(state, listA, N)
    #idea: have a mask for each and keep track of last digit that youve finished.
#shift to that digit
#println("state", state);
    num=0;
    lastNum=-1;
    digitsdone=0;
    for i=1:length(listA)
        temp=0;
        numDigits=listA[i]-lastNum-1;
        #println("numdigits",numDigits, "i, ", listA[i]);
        if(numDigits>0)
            mask=createMask(numDigits, lastNum+1);
            #println("mask",mask);
            #println("lastNum",lastNum);
            temp=(mask&state)>>(lastNum+1);
            #println("temp",temp);
            num|=temp<<digitsdone;
            #println("num",temp);
            digitsdone+=numDigits;
            #println("digits done ",digitsdone);
        end
        lastNum=listA[i];
        println("");
    end
    if(N-1!=listA[length(listA)-1])
        numDigits=N-lastNum-1;
        println("numDigits",numDigits);
        if(numDigits==0)
            return num;
        end
        mask=createMask(numDigits, lastNum+1);
        temp=(mask&state)>>(lastNum+1);
        #println("mask&state",mask&state);
        #println("lastNum+1",lastNum+1);
        num|=temp<<digitsdone;
    end
    #println("final:", num);
    return num;
end




function constructCoefficientMatrix(eigenvector,states,listA, N)
    #listA is a list of the positions of interest
    println("N", N)
    H::Matrix{Float64}=zeros(Float64, 2^length(listA), 2^(N-length(listA)));
    #zeros(Float64, length(list),length(list));
    for i=1:length(states)

        #println("i, ",states[i]);
        #println("listA", listA)
        #1001
        getA=getANum(states[i], listA);
        getB=getBNum(states[i], listA, N);
        #println("A", getA)
        #println("B", getB)

        H[getA+1,getB+1]=eigenvector[i];
    end
    #println(H);
    return H;
end

function getEntanglementEntropy(eigenvector, states, listA, N)
    H=constructCoefficientMatrix(eigenvector, states, listA, N);
    decomp=svd(H);
    tot=0;
    for i=1:length(decomp.S)
        temp=decomp.S[i];
        #println(temp);
        tot+= temp == 0 ? 0 : -(temp^2)*log(temp^2);
    end
    return tot;
end
