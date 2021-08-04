
struct site
    num::Int
    x::Int
    y::Int
end

struct bond
    site1::site
    site2::site
    isNear::Bool;
end


struct refState
    periodicity::Int64;
    state::Int64;
end

struct state
    st::Int64;
    ref::refState;
    shiftsNeeded::Int64;
end

struct refState2d
    state::Int64;
    periodX::Int64;
    periodY::Int64;
    numUniqueSt::Int64;
end

struct state2d
    st::Int64;
    ref::refState2d;
    shiftsXNeeded::Int64;
    shiftsYNeeded::Int64;
end

struct momentum
    px::Int;
    py::Int;
end

struct reflectionRefState
    state::Int64;
    xC::Int64;
    yC::Int64;
    xyC::Int64;
    uniqueReflections::Int64;
end

struct reflectionState
    state::Int64;
    refState::reflectionRefState;
    refXNeeded::Int64;
    refYNeeded::Int64;
end

struct lambda
    lx::Int64;
    ly::Int64;
end

struct square
    num1::site;
    num2::site;
    num3::site;
    num4::site;
end

struct graph
    #(graphNum, numSites, numNearestNeighborBonds, numFarNeighborBonds, numSquares, numPlaquettes, numSubgraphs, latticeConstant, nearBonds, farBonds, squares, plaquettes, siteNumbering)
    num::Int;
    numSites::Int;
    numNearBonds::Int;
    numFarBonds::Int;
    numSquares::Int;
    numPlaquettes::Int;
    numSubgraphs::Int;
    latticeConstant::Float64;
    nearBonds::Vector{bond};
    farBonds::Vector{bond};
    squares::Vector{square};
    plaquettes::Vector{square};
    #1s and 0s labels
    indicies::Vector{Int};
    subgraphList::Vector{Int};

end





function containsSite(i, bond)
    return bond.site1.num==i||bond.site2.num==i;
end
