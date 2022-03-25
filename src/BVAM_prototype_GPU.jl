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
using FiniteDiff

BLAS.set_num_threads(16)
CUDA.allowscalar(false)

#=============================================================
Preconditioners
==============================================================#
function incompletelu(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    Pl = ilu(convert(AbstractMatrix,W), τ = 150.0f0)
  else
    Pl = Plprev
  end
  Pl,nothing
end
Base.eltype(::IncompleteLU.ILUFactorization{Tv,Ti}) where {Tv,Ti} = Tv
#==============================================================
Grid parameters
===============================================================#
dx = 1.0f0
dy = 1.0f0
N = 300
x = range(0, N*dx, length=N)
y = range(0, N*dy, length=N)
tspan = (0.0, Float32(1e3))
#==============================================================
Problem parameters
===============================================================#
println("Setting up problem parameters")
ū = 0.5f0
v̄ = 0.4f0
D = 0.5f0
#D = cu(Float32.(0.1 .+ (0.5-0.1)*exp.(-1/5000*((x.-N*dx/2).^2 .+ (y.-N*dx/2)'.^2))))
#D,_,_,_,_ = init_cond(:ErMnO3_ErOnly, N, M = 9, α=1/80, BC=:Periodic, c01=0.3, dx = dx, f=2/3)
#D = cu(D)
β = -1.0f0 #FIXED
r₁ = 50.0f0
r₂ = 10.0f0
α = -(v̄ - ū*v̄*r₂)./(ū .- ū*v̄^2*r₁)
γ = -α
δ = 1.0f0 #1/(5*pi)^2*(α + β*D)/(2D)

##########################################################
A = [10 1; 1 10] 
#γ_x = kron(A, ones(N÷2,N÷2))
checker = kron(ones(N÷4,N÷4), A)
gauss = 1.0 .+ (10.0-1.0)*exp.(-1/5000*((x[1:N÷2].-N*dx/4).^2 .+ (y[1:N÷2].-N*dx/4)'.^2))
γ_x = [checker 10*ones(N÷2,N÷2); ones(N÷2,N÷2) gauss]
γ_y = [reverse(checker, dims=2) ones(N÷2,N÷2); 10*ones(N÷2,N÷2) (10.0+1.0).-gauss]
#-----------------------------
f = 1e-3
w = 1 
γ_x[N÷2-w+1:N÷2+w,:] .= f
γ_x[:, N÷2-w+1:N÷2+w] .= f
γ_y[N÷2-w+1:N÷2+w,:] .= f
γ_y[:, N÷2-w+1:N÷2+w] .= f
γ_x[1:w, :] .= f
γ_x[end-w+1:end, :] .= f
γ_x[:,1:w] .= f
γ_x[:, end-w+1:end] .= f
γ_y[1:w, :] .= f
γ_y[end-w+1:end, :] .= f
γ_y[:,1:w] .= f
γ_y[:, end-w+1:end] .= f
#----------------
#B = 5 .+ (1-5)*exp.(-1/2500*((x[1:end÷2].-N*dx/4).^2 .+ (y[1:end÷2].-N*dx/4)'.^2))
#A = 5 .+ (9-5)*exp.(-1/2500*((x[1:end÷2].-N*dx/4).^2 .+ (y[1:end÷2].-N*dx/4)'.^2))
#γ_x = cu([A B; B A])
#γ_y = cu(reverse(γ_x, dims=2))
#γ_x = 10+1e-3 .+ 5*sin.(2*pi*5*x) .+ 5*cos.(2*pi*5*y')
#γ_y = 20 .- γ_x
#----------------
γ_x = 2.0f0
γ_y = 1.0f0
#-----------------
A = [10 1; 1 10] 
checker = kron(A, ones(N÷20,N÷20))
γ_x = kron(ones(N÷30, N÷30), checker)
γ_y = reverse(γ_x, dims=2)
#------------------------------------
γ_x = -9(x/2 .+ y'/2)/N.+10
γ_y = reverse(γ_x')
#------------------------------------
γ_x = 5*sin.(10*pi*x.+10*pi*y') .+ (5+1e-3) 
γ_y = 1/2*sin.(10*pi*x.+10*pi*y' .+ pi) .+ (1/2+1e-3) 
##########################################################################################
γᵤ = 1.0f0
γᵥ = 1.0f0
#=============================================================
Initial Conditions
=============================================================#
println("Initial Conditions")
r01 = Float32.(ū .+ rand(Uniform(-ū/100, +ū/100), N, N))
r02 = Float32.(v̄ .+ rand(Uniform(-v̄/100, +v̄/100), N, N))
#r01, x, y, dx, dy = init_cond(:NoisePatches, N, c01 = ū, dx = dx)
#r02, x, y, dx, dy = init_cond(:NoisePatches, N, c01 = v̄, dx = dx)
#r01, x, y, dx, dy = init_cond(:ErMnO3_MnOnly, N, M = 9, α=1/80, BC=:Periodic, c01 =ū, dx = dx)
#r02, x, y, dx, dy = init_cond(:ErMnO3_ErOnly, N, M = 15, α=1, BC=:Neumann, c01 =v̄, dx = dx)
#r01 .= 1 .- r01
r0 = cat(cu(r01), cu(r02), dims=3)
r0_cpu = cat(r01, r02, dims=3)
#=============================================================
Finite differences matrices
==============================================================#
println("Dealing with finite differencing.")
#Q = Neumann0BC(1.0, 1)
Q = PeriodicBC(Float64)
D2 = CenteredDifference(2, 2, 1.0, N)
L2 = D2*Q

A2 = cu(sparse(L2)[1])
A2_cpu = sparse(L2)[1]
#============================================================
Caches
=============================================================#
println("Dealing with caches.")
A2u = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
uA2 = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
D2u = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
A2v = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
vA2 = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
D2v = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
A2u_i = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
u_iA2 = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
u_i = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
A2v_i = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
v_iA2 = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
v_i = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
u_ni = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))
v_ni = DiffEqBase.dualcache(cu(zeros(Float32, N,N)))

cache =(A2u, uA2, D2u, A2v, vA2, D2v, A2u_i, u_iA2, u_i, A2v_i, v_iA2, v_i)
p = (D, δ, γᵤ, γᵥ, α, β, γ, r₁, r₂, cu(γ_x), cu(γ_y), dx, dy, A2, cache)
p_cpu = (D, δ, γᵤ, γᵥ, α, β, γ, r₁, r₂, dx, dy, A2_cpu, cache)
#==========================================================
Advanced versoin of the PDE, for solving
===========================================================#
function advanced_version!(dr,r,p,t)
  D, δ, γᵤ, γᵥ, α, β, γ, r₁, r₂, γ_x, γ_y, dx, dy, A2, cache = p
  A2u, uA2, D2u, A2v, vA2, D2v, A2u_i, u_iA2, u_i, A2v_i, v_iA2, v_i =  cache

  u = @view r[:,:,1]
  v = @view r[:,:,2]
  du = @view dr[:,:,1]
  dv = @view dr[:,:,2]

  A2u = get_tmp(A2u, 0.f0)
  uA2 = get_tmp(uA2, 0.f0)
  D2u = get_tmp(D2u, 0.f0)
  A2v = get_tmp(A2v, 0.f0)
  vA2 = get_tmp(vA2, 0.f0)
  D2v = get_tmp(D2v, 0.f0)
  A2u_i = get_tmp(A2u_i, 0.f0)
  u_iA2 = get_tmp(u_iA2, 0.f0)
  u_i = get_tmp(u_i, 0.f0)
  A2v_i = get_tmp(A2v_i, 0.f0)
  v_iA2 = get_tmp(v_iA2, 0.f0)
  v_i = get_tmp(v_i, 0.f0)

  mul!(A2u,A2,u)
  mul!(uA2,u,cu(A2'))
  mul!(A2v,A2,v)
  mul!(vA2,v,cu(A2'))
  @. u_i = 2.0f0*u*(u-1)*(2.0f0*u-1) - γᵤ*(γ_x*1.0f0/dy^2*A2u + γ_y*1.0f0/dx^2*uA2) # 4 allocations
  @. v_i = 2.0f0*v*(v-1)*(2.0f0*v-1) - γᵥ*(γ_x*1.0f0/dy^2*A2v + γ_y*1.0f0/dx^2*vA2) # 4 allocatoins
  mul!(A2u_i,A2,u_i)
  mul!(u_iA2,u_i,cu(A2'))
  mul!(A2v_i,A2,v_i)
  mul!(v_iA2,v_i,cu(A2'))

  @. D2u = δ*D*(1.0f0/dy^2*A2u_i + 1.0f0/dx^2*u_iA2)
  @. D2v = δ*(1.0f0/dy^2*A2v_i + 1.0f0/dx^2*v_iA2)
  @. du = D2u + α*u*(1.0f0-r₁*v^2) + v*(1.0f0 - r₂*u)
  @. dv = D2v + β*v*(1.0f0+α*r₁/β*u*v) + u*(γ + r₂*v)
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
#@time jac_sparsity = Symbolics.jacobian_sparsity((dr,r)->basic_version!(dr,r,p_cpu,0.f0),dr0,r0_cpu)

#jac_sparsity = load("j500_n0.jld", "jac_sparsity")

#jac = Float32.(sparse(jac_sparsity))
#colors = matrix_colors(jac)
#using FiniteDiff
#sparsecache = FiniteDiff.JacobianCache(r0_cpu, colorvec=colors,sparsity=jac)
#my_basic_version!(dr0, r0) = basic_version!(dr0, r0, p_cpu, 0.f0)
#FiniteDiff.finite_difference_jacobian!(jac, my_basic_version!, r0_cpu, sparsecache, colorvec=colors)
#using Arpack
#@time jac_eig = eigs(jac, nev=2, maxiter=500, ritzvec=false)
#@show jac_eig[1]

#f = ODEFunction(advanced_version!;jac=jac,jac_prototype=Float32.(jac_sparsity))

#==================================================================
Solve the system
==================================================================#
println("Solving")
prob = ODEProblem(advanced_version!,r0,tspan,p)
CUDA.allowscalar(false)
#@time sol = solve(prob, TRBDF2(linsolve=KrylovJL_GMRES(), precs=incompletelu, concrete_jac=true), saveat=range(0.0, stop=tspan[2], length=101), progress=true, progress_steps=1);
@time sol = solve(prob, ROCK2(), saveat=range(0.0, stop=tspan[2], length=101), progress=true, progress_steps=1, reltol=1e-6, abstol=1e-6);

#==================================================================
Plot
==================================================================#
#println("Plotting")
anim = @animate for i in 1:length(sol.t)
  title = "t = $(sol.t[i])"
  #u = kron(ones(3,3), sol.u[i][:,:,1])
  u = Array(sol.u[i])[:,:,1]
  heatmap(x, y, u, c=:grays, aspect_ratio=1, axis=([], false), title=title, colorbar=false)
  #heatmap(u, c=:roma, aspect_ratio=dy/dx, axis=([], false), title=title, colorbar=false)
end

gif(anim, fps=10)