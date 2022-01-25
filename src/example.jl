using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using PatternFormation
using Plots
using Symbolics
using SparseArrays
using SparseDiffTools
using DifferentialEquations, LinearAlgebra
using AlgebraicMultigrid
using LoopVectorization

println("Welcome!")

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
p = [f, k, D₁, D₂, dx, dy, N]

tspan = (0.0, 100.0)
myEquation = GS_Neumann0!

# Jacobian and stuff
du0 = similar(u0)
sparsity_pattern = Symbolics.jacobian_sparsity((du,u)->myEquation(du,u,p,0.0),du0,u0)
jac_sparsity = Float64.(sparse(sparsity_pattern))
colorvec = matrix_colors(jac_sparsity)
ff = ODEFunction(myEquation;jac_prototype=jac_sparsity,colorvec=colorvec)
prob = ODEProblem(ff,u0,tspan,p)

# Preconditioner
function algebraicmultigrid(W,du,u,p,t,newW,Plprev,Prprev,solverdata)
  if newW === nothing || newW
    Pl = aspreconditioner(ruge_stuben(convert(AbstractMatrix,W)))
  else
    Pl = Plprev
  end
  Pl,nothing
end
# Required due to a bug in Krylov.jl: https://github.com/JuliaSmoothOptimizers/Krylov.jl/pull/477
Base.eltype(::AlgebraicMultigrid.Preconditioner) = Float64

# Solve!
println("Solving")
@time sol = solve(prob,KenCarp4(linsolve=KLUFactorization(),precs=algebraicmultigrid), saveat=range(0, stop=tspan[2], length=101), progress=true, progress_steps=1)

# Plot!
anim = @animate for i in 1:length(sol.t)
  tit = "$type-type; f = $f, k = $k; t = $(sol.t[i])" 
  #u = kron(ones(3,3), sol.u[i][:,:,1])
  u = sol.u[i][:,:,1]
  heatmap(u, c=:berlin, aspect_ratio=dx/dy, axis=([], false), title=tit, colorbar=false)
end

gif(anim, fps=10)