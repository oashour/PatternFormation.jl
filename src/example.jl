using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

import LinearAlgebra, OpenBLAS32_jll
LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)

using PatternFormation
using LinearSolve
using Plots
using Symbolics
using SparseArrays
using SparseDiffTools
using OrdinaryDiffEq, LinearAlgebra
using AlgebraicMultigrid
using FLoops
using BenchmarkTools

using IncompleteLU
                      
# Grid and initial conditions
const N = 256
dx = 1/143            
println("Hello")

#u01, x, y, dx, dy = init_cond(:ErMnO3_MnOnly, N, M = 7, α=500, β=3000, BC=:Neumann, r₀=0.001)
#u02, x, y, dx, dy = init_cond(:ErMnO3_ErOnly, N, M = 7, α=500, β=3000, BC=:Neumann, r₀=0.001)
u01, x, y, dx, dy = init_cond(:Square, N, M = 3, α=1000, BC=:Periodic)
u02, x, y, dx, dy = init_cond(:Square, N, M = 3, α=1000, BC=:Periodic)
u01 = 1 .- u01        
u0 = cat(u01, u02, dims=3)
heatmap(u01)
heatmap(u02)

#Parameters
#f = 0.026
#k = 0.051
#type = "β"

f = 0.046
k = 0.065
type = "μ"

D₁ = 2e-5
D₂ = 1e-5
N_threads = 16
BLAS.set_num_threads(N_threads)
ex = ThreadedEx(simd = true)
#ex = ThreadedEx()
p = [f, k, D₁, D₂, dx, dy, N]

tspan = (0.0, 100.0)
myEquation = GS_Periodic!
func(du,u,p,t) = myEquation(du, u, p, t, ex)

# Jacobian and stuff
du0 = similar(u0)
sparsity_pattern = Symbolics.jacobian_sparsity((du,u)->func(du,u,p,0.0),du0,u0)
jac_sparsity = Float64.(sparse(sparsity_pattern))
colorvec = matrix_colors(jac_sparsity)
ff = ODEFunction(func;jac_prototype=jac_sparsity,colorvec=colorvec)
prob = ODEProblem(ff,u0,tspan,p,ex)

function incompletelu(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    Pl = ilu(convert(AbstractMatrix,W), τ = 50)
  else
    Pl = Plprev
  end
  Pl,nothing
end
Base.eltype(::IncompleteLU.ILUFactorization{Tv,Ti}) where {Tv,Ti} = Tv
function algebraicmultigrid(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    Pl = aspreconditioner(ruge_stuben(convert(AbstractMatrix,W)))
  else
    Pl = Plprev
  end
  Pl,nothing
end
function algebraicmultigrid2(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    A = convert(AbstractMatrix,W)
    Pl = AlgebraicMultigrid.aspreconditioner(AlgebraicMultigrid.ruge_stuben(A, presmoother = AlgebraicMultigrid.Jacobi(rand(size(A,1))), postsmoother = AlgebraicMultigrid.Jacobi(rand(size(A,1)))))
  else
    Pl = Plprev
  end
  Pl,nothing
end
Base.eltype(::AlgebraicMultigrid.Preconditioner) = Float64

# Solve!
#println("Solving")
#@time solve(prob,KenCarp4(precs=algebraicmultigrid), saveat=range(0, stop=tspan[2], length=101), progress=true, progress_steps=1)
#println("done")
@benchmark sol = solve(prob, KenCarp47(linsolve=IterativeSolversJL_GMRES()), saveat=range(0, stop=tspan[2], length=101), progress=true, progress_steps=1) samples=1 evals=1 seconds=1
#println("Done!")

# Split
du0 = similar(u0)
func1(du,u,p,t) = GS_Neumann1!(du, u, p, t, ex)
func2(du,u,p,t) = GS_Neumann2!(du, u, p, t, ex)
sparsity_pattern = Symbolics.jacobian_sparsity((du,u)->func1(du,u,p,0.0),du0,u0)
jac_sparsity = Float64.(sparse(sparsity_pattern))
colorvec = matrix_colors(jac_sparsity)
ff = ODEFunction(func1;jac_prototype=jac_sparsity,colorvec=colorvec)
split_prob = SplitODEProblem(ff, GS_Neumann2!, u0, tspan, p)

#println("Solving")
@btime sol = solve(prob,KenCarp47(linsolve=IterativeSolversJL_GMRES(), precs=incompletelu, concrete_jac=true), saveat=range(0, stop=tspan[2], length=101), progress=true, progress_steps=1)


# Plot!
#println("Plotting")
#anim = @animate for i in 1:length(sol.t)
#  tit = "$type-type; f = $f, k = $k; t = $(sol.t[i])" 
#  u = kron(ones(3,3), sol.u[i][:,:,1])
#  #u = sol.u[i][:,:,1]
#  heatmap(u, c=:berlin, aspect_ratio=dx/dy, axis=([], false), title=tit, colorbar=false)
#end

#gif(anim, fps=10)