# For a symmetric real-valued matrix H, returns the givens rotation matrix to zero out the symmetric matrix elements at [i, j] and [j, i]
function calculate_givens_matrix(M, i, j)
    N=size(M, 1)
    # top corner
    a = M[i,i]
    # bottom corner
    c = M[j, j]
    b = M[i, j]

    if(b==0)
        return I
    end

    theta = atan(-2*sqrt((a+sqrt(4*b^2+(a-c)^2)-c)/(sqrt(4*b^2+(a-c)^2))),-(4*b)/(sqrt(4*b^2+(a-c)^2)*sqrt((a+sqrt(4*b^2+(a-c)^2)-c)/(sqrt(4*b^2+(a-c)^2)))))
    Is = copy(sparse(1.0*I, N, N))
    Is[i,i]=cos(theta)
    Is[i, j]=-sin(theta)
    Is[j, i]=sin(theta)
    Is[j,j]=cos(theta)
    return Is
end

