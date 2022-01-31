using MKL

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

#import LinearAlgebra, OpenBLAS32_jll
#LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)

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
N_threads = 1
BLAS.set_num_threads(N_threads)
ex = ThreadedEx(simd = true)
ex = ThreadedEx()
p = [f, k, D₁, D₂, dx, dy, N]

tspan = (0.0, 100.0)
myEquation = GS_Periodic!
func(du,u,p,t) = myEquation(du, u, p, t, ex)

# Solve!
#println("Solving")
#@time sol = solve(prob, ROCK2(), saveat=range(0, stop=tspan[2], length=101), progress=true, progress_steps=1)
#println("Done!")

# Plot!
#println("Plotting")
#anim = @animate for i in 1:length(sol.t)
#  tit = "$type-type; f = $f, k = $k; t = $(sol.t[i])" 
#  u = kron(ones(3,3), sol.u[i][:,:,1])
#  #u = sol.u[i][:,:,1]
#  heatmap(u, c=:berlin, aspect_ratio=dx/dy, axis=([], false), title=tit, colorbar=false)
#end

#gif(anim, fps=10)