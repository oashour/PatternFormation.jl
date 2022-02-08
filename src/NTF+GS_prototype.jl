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
Mᵤ =  1000.0
Mᵥ =  1.0
γᵤ =  1.0
γᵥ =  1.0
ū = 0.73
v̄ = 0.8
F = ū*v̄^2/(1-ū)
k = ū*v̄ - F 
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
r01 = ū .+ rand(Uniform(-ū/100, +ū/100), N, N)
r02 = v̄ .+ rand(Uniform(-v̄/100, +v̄/100), N, N)
x = range(0, N*dx, length=N)
y = range(0, N*dy, length=N)
#r01, x, y, dx, dy = init_cond(:NoisePatches, N, c01 = ū, dx = dx)
#r02, x, y, dx, dy = init_cond(:NoisePatches, N, c01 = v̄, dx = dx)
#r01, x, y, dx, dy = init_cond(:ErMnO3_MnOnly, N, M = 15, α=3000, BC=:Neumann, c01 =ū, f = +0.01, dx = dx)
#r02, x, y, dx, dy = init_cond(:ErMnO3_ErOnly, N, M = 15, α=500, BC=:Neumann, c01 =v̄, f = -0.01, dx = dx)
#r01 .= 1 .- r01
r0 = cat(r01, r02, dims=3)
#=============================================================
Finite differences matrices
==============================================================#
println("Dealing with finite differencing.")
#Q = Neumann0BC(1.0, 1)
Q = PeriodicBC(Float64)
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
p = (Mᵤ, Mᵥ, γᵤ, γᵥ, a, c, dx, dy, A2, cache)
#==========================================================
Advanced versoin of the PDE, for solving
===========================================================#
function advanced_version!(dr,r,p,t)
  Mᵤ, Mᵥ, γᵤ, γᵥ, F, k, dx, dy, A2, cache = p
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

  @. D2u = Mᵤ*(1/dy^2*A2u_i + 1/dx^2*u_iA2)
  @. D2v = Mᵥ*(1/dy^2*A2v_i + 1/dx^2*v_iA2)
  @. du = D2u - u*v^2 .+ F*(1-u)
  @. dv = D2v + u*v^2 .- (F+k)*v
  nothing
end
#=============================================================
Basic Version for computing Jacobian sparsity
==============================================================#
function basic_version!(dr,r,p,t)
  Mᵤ, Mᵥ, γᵤ, γᵥ, F, k, dx, dy, A2, cache = p

  u = @view r[:,:,1]
  v = @view r[:,:,2]
  
  u_n = 2u.*(u.-1).*(2u.-1) .- γᵤ*(1/dy^2*A2*u + 1/dx^2*u*A2')
  v_n = 2v.*(v.-1).*(2v.-1) .- γᵥ*(1/dy^2*A2*v + 1/dx^2*v*A2')
  D2u = Mᵤ*(1/dy^2*A2*u_n + 1/dx^2*u_n*A2')
  D2v = Mᵥ*(1/dy^2*A2*v_n + 1/dx^2*v_n*A2')
  dr[:,:,1] = D2u .- u.*v.^2 .+ F*(1 .- u)
  dr[:,:,2] = D2v .+ u.*v.^2 .- (F+k)*v
  nothing
end
#=================================================================
Compute the Jacobian sparsity pattern
==================================================================#
#println("Jacobian sparsity")
#dr0 = copy(r0)
#@time jac_sparsity = Symbolics.jacobian_sparsity((dr,r)->basic_version!(dr,r,p,0.0),dr0,r0)
f = ODEFunction(advanced_version!;jac_prototype=float.(jac_sparsity))

#==================================================================
Solve the system
==================================================================#
println("Solving")
prob = ODEProblem(f,r0,tspan,p)
@time sol = solve(prob, TRBDF2(linsolve=KrylovJL_GMRES(), precs=incompletelu, concrete_jac=true), saveat=range(0.0, stop=tspan[2], length=101), progress=true, progress_steps=1);


#==================================================================
Plot
==================================================================#
println("Plotting")
anim = @animate for i in 1:length(sol.t)
  title = "t = $(sol.t[i])"
  #u = kron(ones(3,3), sol.u[i][:,:,1])
  u = sol.u[i][:,:,2]
  heatmap(x, y, u, c=:roma, aspect_ratio=1, axis=([], false), title=title, colorbar=false)
end

gif(anim, fps=10)