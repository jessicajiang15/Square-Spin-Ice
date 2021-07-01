
struct site
    num::Int
    x::Int
    y::Int
end

struct bond
    site1::site
    site2::site
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

function containsSite(i, bond)
    return bond.site1.num==i||bond.site2.num==i;
end
