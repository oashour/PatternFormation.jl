using DiffEqOperators
using StaticArrays

N = 10
deriv_order = 2
approx_order = 2
bound_order = 1

Q = Neumann0BC(1.0, bound_order)
D = CenteredDifference(deriv_order, approx_order, 1.0, N)
L = D*Q

conc = sparse(L)[1]


Ax = sparse(Tridiagonal([1.0 for i in 1:N-1],[-2.0 for i in 1:N],[1.0 for i in 1:N-1]))
Ax[2,1] = 2.0
Ax[end-1,end] = 2.0

