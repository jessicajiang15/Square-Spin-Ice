
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

function containsSite(i, bond)
    return bond.site1.num==i||bond.site2.num==i;
end
