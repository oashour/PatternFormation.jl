using MKL

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using OrdinaryDiffEq, LinearAlgebra, Plots, DiffEqOperators, SparseArrays, Distributions
using Symbolics, SparseDiffTools
using PreallocationTools
using PatternFormation
using LinearSolve
using IncompleteLU
using BenchmarkTools
using LazyArrays

BLAS.set_num_threads(16)

#=============================================================
Preconditioners
==============================================================#
function incompletelu(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    Pl = ilu(convert(AbstractMatrix,W), τ = 150.0)
  else
    Pl = Plprev
  end
  Pl,nothing
end
Base.eltype(::IncompleteLU.ILUFactorization{Tv,Ti}) where {Tv,Ti} = Tv
#==============================================================
Problem parameters
===============================================================#
println("Setting up problem parameters")
ū = 0.88
v̄ = 0.37
D = 0.5
β = -1.0 #FIXEf0
r₁ = 3.0
r₂ = 10.0
α = -(v̄ - ū*v̄*r₂)/(ū - ū*v̄^2*r₁)
γ = -α
δ = 1.0; #1/(5*pi)^2*(α + β*D)/(2D)
γᵤ =  10.0
γᵥ =   1.0
#==============================================================
Grid parameters
===============================================================#
dx = 1
dy = 1
N = 200
tspan = (0.0, 1e3)
#=============================================================
Initial Conditions
=============================================================#
println("Initial Conditions")
#r01 = rand(Uniform(-0.5, +0.5), N, N)
#r02 = rand(Uniform(-0.5, +0.5), N, N)
r01 = ū .+ rand(Uniform(-ū/100, +ū/100), N, N)
r02 = v̄ .+ rand(Uniform(-v̄/100, +v̄/100), N, N)
x = range(0, N*dx, length=N)
y = range(0, N*dy, length=N)
#r01, x, y, dx, dy = init_cond(:NoisePatches, N, c01 = ū, dx = dx)
#r02, x, y, dx, dy = init_cond(:NoisePatches, N, c01 = v̄, dx = dx)
#r01, x, y, dx, dy = init_cond(:ErMnO3_MnOnly, N, M = 15, α=0.3*4, BC=:Periodic, c01 =ū, f = 0.01, dx = dx)
#r02, x, y, dx, dy = init_cond(:ErMnO3_ErOnly, N, M = 15, α=0.06*4, BC=:Periodic, c01 =v̄, f = 0.01, dx = dx)
#r01 .= 1 .- r01
r0 = cat(r01, r02, dims=3)
#=============================================================
Finite differences matrices
==============================================================#
println("Dealing with finite differencing.")
Q = Neumann0BC(1.0, 1)
#Q = PeriodicBC(Float64)
D2 = CenteredDifference(2, 2, 1.0, N)
L2 = D2*Q

A2 = sparse(L2)[1]
#============================================================
Caches
=============================================================#
println("Dealing with caches.")
A2u = DiffEqBase.dualcache(zeros(N,N))
uA2 = DiffEqBase.dualcache(zeros(N,N))
D2u = DiffEqBase.dualcache(zeros(N,N))
A2v = DiffEqBase.dualcache(zeros(N,N))
vA2 = DiffEqBase.dualcache(zeros(N,N))
D2v = DiffEqBase.dualcache(zeros(N,N))
A2u_i = DiffEqBase.dualcache(zeros(N,N))
u_iA2 = DiffEqBase.dualcache(zeros(N,N))
u_i = DiffEqBase.dualcache(zeros(N,N))
A2v_i = DiffEqBase.dualcache(zeros(N,N))
v_iA2 = DiffEqBase.dualcache(zeros(N,N))
v_i = DiffEqBase.dualcache(zeros(N,N))
u_ni = DiffEqBase.dualcache(zeros(N,N))
v_ni = DiffEqBase.dualcache(zeros(N,N))

cache =(A2u, uA2, D2u, A2v, vA2, D2v, A2u_i, u_iA2, u_i, A2v_i, v_iA2, v_i)
p = (D, δ, γᵤ, γᵥ, α, β, γ, r₁, r₂, dx, dy, A2, cache)
#==========================================================
Advanced versoin of the PDE, for solving
===========================================================#
function advanced_version!(dr,r,p,t)
  D, δ, γᵤ, γᵥ, α, β, γ, r₁, r₂, dx, dy, A2, cache = p
  A2u, uA2, D2u, A2v, vA2, D2v, A2u_i, u_iA2, u_i, A2v_i, v_iA2, v_i =  cache

  u = @view r[:,:,1]
  v = @view r[:,:,2]
  du = @view dr[:,:,1]
  dv = @view dr[:,:,2]

  A2u = get_tmp(A2u, first(r)*t)
  uA2 = get_tmp(uA2, first(r)*t)
  D2u = get_tmp(D2u, first(r)*t)
  A2v = get_tmp(A2v, first(r)*t)
  vA2 = get_tmp(vA2, first(r)*t)
  D2v = get_tmp(D2v, first(r)*t)
  A2u_i = get_tmp(A2u_i, first(r)*t)
  u_iA2 = get_tmp(u_iA2, first(r)*t)
  u_i = get_tmp(u_i, first(r)*t)
  A2v_i = get_tmp(A2v_i, first(r)*t)
  v_iA2 = get_tmp(v_iA2, first(r)*t)
  v_i = get_tmp(v_i, first(r)*t)
  
  mul!(A2u,A2,u)
  mul!(uA2,u,A2')
  mul!(A2v,A2,v)
  mul!(vA2,v,A2')
  @. u_i = 2u*(u-1)*(2u-1) - γᵤ*(1/dy^2*A2u + 1/dx^2*uA2) # 4 allocations
  @. v_i = 2v*(v-1)*(2v-1) - γᵥ*(1/dy^2*A2v + 1/dx^2*vA2) # 4 allocatoins
  mul!(A2u_i,A2,u_i)
  mul!(u_iA2,u_i,A2')
  mul!(A2v_i,A2,v_i)
  mul!(v_iA2,v_i,A2')

  @. D2u = δ*D*(1/dy^2*A2u_i + 1/dx^2*u_iA2)
  @. D2v = δ*(1/dy^2*A2v_i + 1/dx^2*v_iA2)
  #@. D2u = δ*D*(1/dy^2*A2u + 1/dx^2*uA2)
  #@. D2v = δ*(1/dy^2*A2v + 1/dx^2*vA2)
  @. du = D2u + α*u*(1-r₁*v^2) + v*(1 - r₂*u)
  @. dv = D2v + β*v*(1+α*r₁/β*u*v) + u*(γ + r₂*v)
  nothing
end
#=============================================================
Basic Version for computing Jacobian sparsity
==============================================================#
function basic_version!(dr,r,p,t)
  D, δ, γᵤ, γᵥ, α, β, γ, r₁, r₂, dx, dy, A2, cache = p

  u = @view r[:,:,1]
  v = @view r[:,:,2]
  
  u_n = 2u.*(u.-1).*(2u.-1) .- γᵤ*(1/dy^2*A2*u + 1/dx^2*u*A2')
  v_n = 2v.*(v.-1).*(2v.-1) .- γᵥ*(1/dy^2*A2*v + 1/dx^2*v*A2')
  D2u = D*δ*(1/dy^2*A2*u_n + 1/dx^2*u_n*A2')
  D2v = δ*(1/dy^2*A2*v_n + 1/dx^2*v_n*A2')
  #D2u = D*δ*(1/dy^2*A2*u + 1/dx^2*u*A2')
  #D2v = δ*(1/dy^2*A2*v + 1/dx^2*v*A2')
  dr[:,:,1] = D2u .+ α*u.*(1 .- r₁*v.^2) .+ v.*(1 .- r₂*u)
  dr[:,:,2] = D2v .+ β*v.*(1 .+ α*r₁/β*u.*v) .+ u.*(γ .+ r₂*v)
  nothing
end
#=================================================================
Compute the Jacobian sparsity pattern
==================================================================#
println("Jacobian sparsity")
dr0 = copy(r0)
@time jac_sparsity = Symbolics.jacobian_sparsity((dr,r)->basic_version!(dr,r,p,0.0),dr0,r0)
f = ODEFunction(advanced_version!;jac_prototype=float.(jac_sparsity))

#==================================================================
Solve the system
==================================================================#
println("Solving")
prob = ODEProblem(f,r0,tspan,p)
@time sol = solve(prob, TRBDF2(linsolve=KrylovJL_GMRES(), precs=incompletelu, concrete_jac=true), saveat=range(0.0, stop=tspan[2], length=101), progress=true, progress_steps=1);
#@time sol = solve(prob, ROCK2(), saveat=range(0.0, stop=tspan[2], length=101), progress=true, progress_steps=1);


#==================================================================
Plot
==================================================================#
println("Plotting")
anim = @animate for i in 1:length(sol.t)
  title = "t = $(sol.t[i])"
  #u = kron(ones(3,3), sol.u[i][:,:,1])
  u = sol.u[i][:,:,2]
  heatmap(x, y, u, c=:roma, aspect_ratio=1, axis=([], false), title=title, colorbar=false)
  #heatmap(u, c=:roma, aspect_ratio=dy/dx, axis=([], false), title=title, colorbar=false)
end

gif(anim, fps=10)