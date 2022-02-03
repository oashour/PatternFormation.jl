using MKL

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

#import LinearAlgebra, OpenBLAS32_jll
#LinearAlgebra.BLAS.lbt_forward(OpenBLAS32_jll.libopenblas_path)

using PatternFormation
using LinearSolve
using Plots
using OrdinaryDiffEq, LinearAlgebra
using FLoops
using BenchmarkTools
                      
# Grid and initial conditions
const N = 2000
myEquation = GS_Neumann0!
tspan = (0.0, 10000.0)

f = 0.090
k = 0.057
type = "σ"

D₁ = 2e-5
D₂ = 1e-5
N_threads = 1
BLAS.set_num_threads(N_threads)
ex = ThreadedEx(simd = true)

println("Hello")

u01, x, y, dx, dy = init_cond(:ErMnO3_MnOnly, N, M = 21, α=3000, γ = 500, BC=:Neumann, r₀=0.01, n_vd = 100, n_id = 100)
u02, x, y, dx, dy = init_cond(:ErMnO3_ErOnly, N, M = 21, α=500, γ = 3000, BC=:Neumann, r₀=0.01, n_vd = 100, n_id = 100)
#u01, x, y, dx, dy = init_cond(:ErMnO3_MnOnly, N, M = 21, α=3000, BC=:Neumann, r₀=0.01)
#u02, x, y, dx, dy = init_cond(:ErMnO3_ErOnly, N, M = 21, α=500, BC=:Neumann, r₀=0.01)
#u01, x, y, dx, dy = init_cond(:ErMnO3_MnOnly, N, M = 21, α=3000, BC=:Periodic)
#u02, x, y, dx, dy = init_cond(:ErMnO3_ErOnly, N, M = 21, α=500, BC=:Periodic)
u01 = 1 .- u01        
u0 = cat(u01, u02, dims=3)
#heatmap(x, y, u01')
#heatmap(x, y, u02')
p = [f, k, D₁, D₂, dx, dy, N]

func(du,u,p,t) = myEquation(du, u, p, t, ex)
prob = ODEProblem(func, u0, tspan, p)

# Solve!
println("Solving")
@time sol = solve(prob, ROCK2(), saveat=range(0, stop=tspan[2], length=101), progress=true, progress_steps=1)
println("Done!")

# Plot!
println("Plotting")
anim = @animate for i in 1:length(sol.t)
  tit = "$type-type; f = $f, k = $k; t = $(sol.t[i])" 
  #u = kron(ones(3,3), sol.u[i][:,:,1])
  u = sol.u[i][:,:,1]
  heatmap(u, c=:grays, aspect_ratio=dx/dy, axis=([], false), title=tit, colorbar=false, clims=(minimum(u),maximum(u)))
end

gif(anim, fps=10)