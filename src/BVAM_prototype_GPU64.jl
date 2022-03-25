using MKL

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using OrdinaryDiffEq, LinearAlgebra, Plots, DiffEqOperators, SparseArrays, Distributions, CUDA
using Symbolics, SparseDiffTools
using PreallocationTools
using PatternFormation
using LinearSolve
using IncompleteLU
using BenchmarkTools
using LazyArrays

BLAS.set_num_threads(16)
CUDA.allowscalar(false)

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
β = -1.0 #FIXE
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
dx = Float64(1/5)
dy = Float64(1/5)
N = 500
tspan = (0.0, Float64(1e2))
#=============================================================
Initial Conditions
=============================================================#
println("Initial Conditions")
r01 = Float64.(ū .+ rand(Uniform(-ū/100, +ū/100), N, N))
r02 = Float64.(v̄ .+ rand(Uniform(-v̄/100, +v̄/100), N, N))
x = range(0, N*dx, length=N)
y = range(0, N*dy, length=N)
#r01, x, y, dx, dy = init_cond(:NoisePatches, N, c01 = ū, dx = dx)
#r02, x, y, dx, dy = init_cond(:NoisePatches, N, c01 = v̄, dx = dx)
#r01, x, y, dx, dy = init_cond(:ErMnO3_MnOnly, N, M = 15, α=3000, BC=:Neumann, c01 =ū, f = +0.01, dx = dx)
#r02, x, y, dx, dy = init_cond(:ErMnO3_ErOnly, N, M = 15, α=500, BC=:Neumann, c01 =v̄, f = -0.01, dx = dx)
#r01 .= 1 .- r01
r0 = cat(CuArray(r01), CuArray(r02), dims=3)
r0_cpu = cat(r01, r02, dims=3)
#=============================================================
Finite differences matrices
==============================================================#
println("Dealing with finite differencing.")
Q = Neumann0BC(1.0, 1)
#Q = PeriodicBC(Float64)
D2 = CenteredDifference(2, 2, 1.0, N)
L2 = D2*Q

A2 = CuArray(sparse(L2)[1])
A2_cpu = sparse(L2)[1]
#============================================================
Caches
=============================================================#
println("Dealing with caches.")
A2u = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
uA2 = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
D2u = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
A2v = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
vA2 = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
D2v = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
A2u_i = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
u_iA2 = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
u_i = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
A2v_i = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
v_iA2 = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
v_i = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
u_ni = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))
v_ni = DiffEqBase.dualcache(CuArray(zeros(Float64, N,N)))

cache =(A2u, uA2, D2u, A2v, vA2, D2v, A2u_i, u_iA2, u_i, A2v_i, v_iA2, v_i)
p = (D, δ, γᵤ, γᵥ, α, β, γ, r₁, r₂, dx, dy, A2, cache)
p_cpu = (D, δ, γᵤ, γᵥ, α, β, γ, r₁, r₂, dx, dy, A2_cpu, cache)
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

  A2u = get_tmp(A2u, 0.0)
  uA2 = get_tmp(uA2, 0.0)
  D2u = get_tmp(D2u, 0.0)
  A2v = get_tmp(A2v, 0.0)
  vA2 = get_tmp(vA2, 0.0)
  D2v = get_tmp(D2v, 0.0)
  A2u_i = get_tmp(A2u_i, 0.0)
  u_iA2 = get_tmp(u_iA2, 0.0)
  u_i = get_tmp(u_i, 0.0)
  A2v_i = get_tmp(A2v_i, 0.0)
  v_iA2 = get_tmp(v_iA2, 0.0)
  v_i = get_tmp(v_i, 0.0)

  mul!(A2u,A2,u)
  mul!(uA2,u,CuArray(A2'))
  mul!(A2v,A2,v)
  mul!(vA2,v,CuArray(A2'))
  @. u_i = 2.0*u*(u-1)*(2.0*u-1) - γᵤ*(1.0/dy^2*A2u + 1.0/dx^2*uA2) # 4 allocations
  @. v_i = 2.0*v*(v-1)*(2.0*v-1) - γᵥ*(1.0/dy^2*A2v + 1.0/dx^2*vA2) # 4 allocatoins
  mul!(A2u_i,A2,u_i)
  mul!(u_iA2,u_i,CuArray(A2'))
  mul!(A2v_i,A2,v_i)
  mul!(v_iA2,v_i,CuArray(A2'))

  @. D2u = δ*D*(1.0/dy^2*A2u_i + 1.0/dx^2*u_iA2)
  @. D2v = δ*(1.0/dy^2*A2v_i + 1.0/dx^2*v_iA2)
  @. du = D2u + α*u*(1.0-r₁*v^2) + v*(1.0 - r₂*u)
  @. dv = D2v + β*v*(1.0+α*r₁/β*u*v) + u*(γ + r₂*v)
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
  dr[:,:,1] = D2u .+ α*u.*(1 .- r₁*v.^2) .+ v.*(1 .- r₂*u)
  dr[:,:,2] = D2v .+ β*v.*(1 .+ α*r₁/β*u.*v) .+ u.*(γ .+ r₂*v)
  nothing
end
#=================================================================
Compute the Jacobian sparsity pattern
==================================================================#
#println("Jacobian sparsity")
#dr0 = copy(r0_cpu)
#@time jac_sparsity = Symbolics.jacobian_sparsity((dr,r)->basic_version!(dr,r,p_cpu,0.),dr0,r0_cpu)
#f = ODEFunction(advanced_version!;jac_prototype=Float64.(jac_sparsity))

#==================================================================
Solve the system
==================================================================#
println("Solving")
prob = ODEProblem(advanced_version!,r0,tspan,p)
CUDA.allowscalar(false)
#@time sol = solve(prob, TRBDF2(linsolve=KrylovJL_GMRES(), precs=incompletelu, concrete_jac=true), saveat=range(0.0, stop=tspan[2], length=101), progress=true, progress_steps=1);
@time sol = solve(prob, ROCK2(), saveat=range(0.0, stop=tspan[2], length=101), progress=true, progress_steps=1);
#@time sol = solve(prob, ROCK2(), saveat=[tspan[2]], progress=true, progress_steps=1);

#==================================================================
Plot
==================================================================#
#println("Plotting")
anim = @animate for i in 1:length(sol.t)
  title = "t = $(sol.t[i])"
  #u = kron(ones(3,3), sol.u[i][:,:,1])
  u = Array(sol.u[i])[:,:,2]
  heatmap(x, y, u, c=:roma, aspect_ratio=1, axis=([], false), title=title, colorbar=false)
  #heatmap(u, c=:roma, aspect_ratio=dy/dx, axis=([], false), title=title, colorbar=false)
end

gif(anim, fps=10)
#title = "t = $(tspan[2])"
#u = Array(sol.u[end])[:,:,2]
#heatmap(u, c=:roma, aspect_ratio=dy/dx, axis=([], false), title=title, colorbar=true)