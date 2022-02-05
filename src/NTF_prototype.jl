using MKL

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using OrdinaryDiffEq, LinearAlgebra, Plots, DiffEqOperators, SparseArrays, Distributions
using Symbolics, SparseDiffTools
using PreallocationTools
using PatternFormation

using LinearSolve

BLAS.set_num_threads(1)

# Generate the constants
println("Setting up problem parameters")
Mᵤ = 5*1e-5
Mᵥ = 1e-5
γᵤ = 1e-5
γᵥ = 1e-5
a = 0.062 # F
c = 0.061 # k
type = "π"

dx = 1/50
dy = 1/50 # Can't be too small? Instability Detected
N = 1000 
tspan = (0.0, 1000.0)

# IC
println("Initial Conditions")
r01 = 0.01 .+ rand(Uniform(-0.01/100, +0.01/100), N, N)
r02 = 0.25 .+ rand(Uniform(-0.25/100, +0.25/100), N, N)
#r01, x, y, dx, dy = init_cond(:NoisePatches, N, c01 = 0.01, dx = dx)
#r02, x, y, dx, dy = init_cond(:NoisePatches, N, c01 = 0.25, dx = dx)
#r01, x, y, dx, dy = init_cond(:ErMnO3_MnOnly, N, M = 9, α=3000, BC=:Neumann)
#r02, x, y, dx, dy = init_cond(:ErMnO3_ErOnly, N, M = 9, α=500, BC=:Neumann)
r01 .= 1 .- r01
r0 = cat(r01, r02, dims=3)

# Generate matrices and stuff
println("Dealing with finite differencing.")
# See https://www.ljll.math.upmc.fr/frey/ftp/finite-differences.pdf to compare to the general scheme
#Q = Neumann0BC(1.0, 1)
Q = PeriodicBC(Float64)
D4 = CenteredDifference(4, 2, 1.0, N)
D2 = CenteredDifference(2, 2, 1.0, N)
D1 = CenteredDifference(1, 2, 1.0, N)
L4 = D4*Q
L2 = D2*Q
L1 = D1*Q

A4 = sparse(L4)[1]
A2 = sparse(L2)[1]
A1 = sparse(L1)[1]

#A2 = sparse(Tridiagonal([1.0 for i in 1:N-1],[-2.0 for i in 1:N],[1.0 for i in 1:N-1]))
#A2[2,1] = 2.0
#A2[end-1,end] = 2.0

println("Dealing with caches.")
A4u = DiffEqBase.dualcache(zeros(N,N))
uA4 = DiffEqBase.dualcache(zeros(N,N))
D4u = DiffEqBase.dualcache(zeros(N,N))
A4v = DiffEqBase.dualcache(zeros(N,N))
vA4 = DiffEqBase.dualcache(zeros(N,N))
D4v = DiffEqBase.dualcache(zeros(N,N))
A2u = DiffEqBase.dualcache(zeros(N,N))
uA2 = DiffEqBase.dualcache(zeros(N,N))
D2u = DiffEqBase.dualcache(zeros(N,N))
A2v = DiffEqBase.dualcache(zeros(N,N))
vA2 = DiffEqBase.dualcache(zeros(N,N))
D2v = DiffEqBase.dualcache(zeros(N,N))
A1u = DiffEqBase.dualcache(zeros(N,N))
uA1 = DiffEqBase.dualcache(zeros(N,N))
D1u = DiffEqBase.dualcache(zeros(N,N))
A1v = DiffEqBase.dualcache(zeros(N,N))
vA1 = DiffEqBase.dualcache(zeros(N,N))
D1v = DiffEqBase.dualcache(zeros(N,N))

cache =  [A4u, uA4, D4u, A4v, vA4, D4v, A2u, uA2, D2u, A2v, vA2, D2v, A1u, uA1, D1u, A1v, vA1, D1v]
p = [Mᵤ, Mᵥ, γᵤ, γᵥ, a, c, dx, dy, A1, A2, A4, cache]

#function basic_version_GS!(dr,r,p,t)
#  Mᵤ, Mᵥ, γᵤ, γᵥ, a, c, dx, dy, A1, A2, A4 = p
#
#  u = @view r[:,:,1]
#  v = @view r[:,:,2]
#  D2u = Mᵤ*(1/dy^2*A2*u + 1/dx^2*u*A2')
#  D2v = Mᵥ*(1/dy^2*A2*v + 1/dx^2*v*A2')
#  dr[:,:,1] = D2u .- u.*v.^2 + a*(1 .- u)
#  dr[:,:,2] = D2v .+ u.*v.^2 - (a+c)*v
#end

function advanced_version!(dr,r,p,t)
  Mᵤ, Mᵥ, γᵤ, γᵥ, F, k, dx, dy, A1, A2, A4, cache = p
  A4u, uA4, D4u, A4v, vA4, D4v, A2u, uA2, D2u, A2v, vA2, D2v, A1u, uA1, D1u, A1v, vA1, D1v = cache

  u = @view r[:,:,1]
  v = @view r[:,:,2]
  du = @view dr[:,:,1]
  dv = @view dr[:,:,2]

  A4u = get_tmp(A4u, first(r)*t)
  uA4 = get_tmp(uA4, first(r)*t)
  D4u = get_tmp(D4u, first(r)*t)
  A4v = get_tmp(A4v, first(r)*t)
  vA4 = get_tmp(vA4, first(r)*t)
  D4v = get_tmp(D4v, first(r)*t)
  A2u = get_tmp(A2u, first(r)*t)
  uA2 = get_tmp(uA2, first(r)*t)
  D2u = get_tmp(D2u, first(r)*t)
  A2v = get_tmp(A2v, first(r)*t)
  vA2 = get_tmp(vA2, first(r)*t)
  D2v = get_tmp(D2v, first(r)*t)
  A1u = get_tmp(A1u, first(r)*t)
  uA1 = get_tmp(uA1, first(r)*t)
  D1u = get_tmp(D1u, first(r)*t) 
  A1v = get_tmp(A1v, first(r)*t)
  vA1 = get_tmp(vA1, first(r)*t)
  D1v = get_tmp(D1v, first(r)*t)
  
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
  @. D1u = 12Mᵤ*(2u .- 1).*(1/dy*A1u + 1/dx*uA1).^2
  @. D1v = 12Mᵥ*(2v .- 1).*(1/dy*A1v + 1/dx*vA1).^2
  @. du = D4u + D2u + D1u + c - a*u*v*v + F*(1-u)
  @. dv = D4v + D2v + D1v + c + a*u*v*v - (F+k)*v
  nothing
end

function basic_version!(dr,r,p,t)
  Mᵤ, Mᵥ, γᵤ, γᵥ, F, k, dx, dy, A1, A2, A4 = p

  u = @view r[:,:,1]
  v = @view r[:,:,2]
  D4u = -γᵤ*Mᵤ*(1/dy^4*A4*u + 1/dx^4*u*A4')
  D4v = -γᵥ*Mᵥ*(1/dy^4*A4*v + 1/dx^4*v*A4')
  D2u = 2Mᵤ*(6u.^2 .- 6u .+ 1).*(1/dy^2*A2*u + 1/dx^2*u*A2')
  D2v = 2Mᵥ*(6v.^2 .- 6v .+ 1).*(1/dy^2*A2*v + 1/dx^2*v*A2')
  D1u = 12Mᵤ*(2u .- 1).*(1/dy*A1*u + 1/dx*u*A1').^2
  D1v = 12Mᵥ*(2v .- 1).*(1/dy*A1*v + 1/dx*v*A1').^2
  dr[:,:,1] = D4u .+ D2u .+ D1u .- u.*v.^2 + F*(1 .- u)
  dr[:,:,2] = D4v .+ D2v .+ D1v .+ u.*v.^2 - (F+k)*v
  nothing
end

# Jacobian stuff
println("Jacobian sparsity")
#dr0 = copy(r0)
#@time jac_sparsity = Symbolics.jacobian_sparsity((dr,r)->basic_version!(dr,r,p,0.0),dr0,r0)
f = ODEFunction(advanced_version!;jac_prototype=float.(jac_sparsity))

println("Solving")
prob = ODEProblem(f,r0,tspan,p)
@time sol = solve(prob, ROCK2(), saveat=range(0.0, stop=tspan[2], length=101), progress=true, progress_steps=1)

println("Plotting")
anim = @animate for i in 1:length(sol.t)
  title = "$type-type; f = $a, k = $c; t = $(sol.t[i])" 
  #u = kron(ones(3,3), sol.u[i][:,:,1])
  u = sol.u[i][:,:,2]
  heatmap(u, c=:berlin, aspect_ratio=dx/dy, axis=([], false), title=title, colorbar=false)
end

gif(anim, fps=10)