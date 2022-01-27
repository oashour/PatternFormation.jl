using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using LinearSolve
using PatternFormation
using Plots
using Symbolics
using SparseArrays
using SparseDiffTools
using OrdinaryDiffEq, LinearAlgebra
using AlgebraicMultigrid
using FLoops
using IterativeSolvers
using BenchmarkTools
                      
# Grid and initial conditions
const N = 256         
dx = 1/143            
u01, x, y, dx, dy = init_cond(:ErMnO3_MnOnly, N, M = 7, α=500, β=3000, BC=:Neumann, r₀=0.001)
u02, x, y, dx, dy = init_cond(:ErMnO3_ErOnly, N, M = 7, α=500, β=3000, BC=:Neumann, r₀=0.001)
u01 = 1 .- u01        
u0 = cat(u01, u02, dims=3)

#Parameters
f = 0.046
k = 0.065
type = "μ"
D₁ = 2e-5
D₂ = 1e-5
N_threads = 1
BLAS.set_num_threads(N_threads)
#ex = ThreadedEx(simd = Val(true))
ex = ThreadedEx()
p = [f, k, D₁, D₂, dx, dy, N]

tspan = (0.0, 1.0)
myEquation = GS_Neumann0!
func(du,u,p,t) = myEquation(du, u, p, t, ex)

# Jacobian and stuff
du0 = similar(u0)
sparsity_pattern = Symbolics.jacobian_sparsity((du,u)->func(du,u,p,0.0),du0,u0)
jac_sparsity = Float64.(sparse(sparsity_pattern))
colorvec = matrix_colors(jac_sparsity)
ff = ODEFunction(func;jac_prototype=jac_sparsity,colorvec=colorvec)
prob = ODEProblem(ff,u0,tspan,p)

using IncompleteLU
function incompletelu(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    Pl = ilu(convert(AbstractMatrix,W), τ = 50.0)
  else
    Pl = Plprev
  end
  Pl,nothing
end
# Required due to a bug in Krylov.jl: https://github.com/JuliaSmoothOptimizers/Krylov.jl/pull/477
Base.eltype(::IncompleteLU.ILUFactorization{Tv,Ti}) where {Tv,Ti} = Tv

# Solve!
println("Solving")
#@profview for i in 1:100 solve(prob,KenCarp4(linsolve=linsolve=KrylovJL_GMRES(), precs=incompletelu, concrete_jac=true), saveat=range(0, stop=tspan[2], length=101)) end
#println("done")
@time sol = solve(prob,KenCarp47(linsolve=KrylovJL_GMRES(), precs=incompletelu, concrete_jac=true), saveat=range(0, stop=tspan[2], length=101), progress=true, progress_steps=1)

# Split
#split_prob = SplitODEProblem(GS_Neumann1!, GS_Neumann2!, u0, tspan, p)

#println("Solving")
#@time sol = solve(prob,KenCarp4(), saveat=range(0, stop=tspan[2], length=101))


# Plot!
#anim = @animate for i in 1:length(sol.t)
#  tit = "$type-type; f = $f, k = $k; t = $(sol.t[i])" 
#  #u = kron(ones(3,3), sol.u[i][:,:,1])
#  u = sol.u[i][:,:,1]
#  heatmap(u, c=:berlin, aspect_ratio=dx/dy, axis=([], false), title=tit, colorbar=false)
#end
#
#gif(anim, fps=10)