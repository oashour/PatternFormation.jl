using OrdinaryDiffEq, LinearAlgebra, Plots, DiffEqOperators, SparseArrays, Distributions
using Symbolics, SparseDiffTools


# Generate the constants
Mᵤ = 1e-5
Mᵥ = 1e-5
γᵤ = 1e-3
γᵥ = 1e-3
a = 0.5
c = 0.00125
dx = 1/100
dy = 1/100

N = 1000

# See https://www.ljll.math.upmc.fr/frey/ftp/finite-differences.pdf to compare to the general scheme
Q = Neumann0BC(1.0, 1)
#Q = PeriodicBC(Float64)
D4 = CenteredDifference(4, 2, 1.0, N)
D2 = CenteredDifference(2, 2, 1.0, N)
D1 = CenteredDifference(1, 2, 1.0, N)
L4 = D4*Q
L2 = D2*Q
L1 = D1*Q

A4 = sparse(L4)[1]
A2 = sparse(L2)[1]
A1 = sparse(L1)[1]


A4u = zeros(N,N)
uA4 = zeros(N,N)
D4u = zeros(N,N)
A4v = zeros(N,N)
vA4 = zeros(N,N)
D4v = zeros(N,N)
A2u = zeros(N,N)
uA2 = zeros(N,N)
D2u = zeros(N,N)
A2v = zeros(N,N)
vA2 = zeros(N,N)
D2v = zeros(N,N)
A1u = zeros(N,N)
uA1 = zeros(N,N)
D1u = zeros(N,N)
A1v = zeros(N,N)
vA1 = zeros(N,N)
D1v = zeros(N,N)

cache =  [A4u, uA4, D4u, A4v, vA4, D4v, A2u, uA2, D2u, A2v, vA2, D2v, A1u, uA1, D1u, A1v, vA1, D1v]
p = [Mᵤ, Mᵥ, γᵤ, γᵥ, a, c, dx, dy, A1, A2, A4, cache]

function basic_version!(dr,r,p,t)
  Mᵤ, Mᵥ, γᵤ, γᵥ, a, c, dx, dy, A1, A2, A4, cache = p
  A4u, uA4, D4u, A4v, vA4, D4v, A2u, uA2, D2u, A2v, vA2, D2v, A1u, uA1, D1u, A1v, vA1, D1v = cache

  u = @view r[:,:,1]
  v = @view r[:,:,2]
  du = @view dr[:,:,1]
  dv = @view dr[:,:,2]

  mul!(A4u,A4,u)
  mul!(uA4,u,A4')
  mul!(A4v,A4,v)
  mul!(vA4,v,A4')
  mul!(A2u,A2,u)
  mul!(uA2,u,A2')
  mul!(A2v,A2,v)
  mul!(vA2,v,A2')
  mul!(A1u,A1,u)
  mul!(uA1,u,A1')
  mul!(A1v,A1,v)
  mul!(vA1,v,A1')

  @. D4u = -γᵤ*Mᵤ*(1/dy^4*A4u + 1/dx^4*uA4)
  @. D4v = -γᵥ*Mᵥ*(1/dy^4*A4v + 1/dx^4*vA4)
  @. D2u = 2Mᵤ*(6u.^2 .- 6u .+ 1).*(1/dy^2*A2u + 1/dx^2*uA2)
  @. D2v = 2Mᵥ*(6v.^2 .- 6v .+ 1).*(1/dy^2*A2v + 1/dx^2*vA2)
  @. D1u = 12Mᵤ*(2u .- 1).*(1/dy*A1u + 1/dx*uA1)
  @. D1v = 12Mᵥ*(2v .- 1).*(1/dy*A1v + 1/dx*vA1)
  @. du = D4u + D2u + D1u + c - a*v*u 
  @. dv = D4v + D2v + D1v + c - a*v*u
end


r01 = 0.01 .+ rand(Uniform(-0.01/100, +0.01/100), N, N)
r02 = 0.25 .+ rand(Uniform(-0.25/100, +0.25/100), N, N)
r0 = cat(r01, r02, dims=3)

#dr0 = copy(r0)
#jac_sparsity = Symbolics.jacobian_sparsity((dr,r)->basic_version!(dr,r,p,0.0),dr0,r0)
#f = ODEFunction(basic_version!;jac_prototype=float.(jac_sparsity))

tspan = (0.0, 5000.0)
prob = ODEProblem(basic_version!,r0,tspan,p)
@time sol = solve(prob, ROCK2(), saveat=range(0.0, stop=tspan[2], length=101), progress=true, progress_steps=1)

println("Plotting")
anim = @animate for i in 1:length(sol.t)
  tit = "t = $(sol.t[i])" 
  #u = kron(ones(3,3), sol.u[i][:,:,1])
  u = sol.u[i][:,:,1]
  heatmap(u, c=:grays, aspect_ratio=dx/dy, axis=([], false), title=tit, colorbar=false, clims=(minimum(u),maximum(u)))
end

gif(anim, fps=10)