using DiffEqOperators, Setfield, Parameters
using LinearAlgebra, Plots, SparseArrays
using BifurcationKit

plotsol(x, Nx=Nx, Ny=Ny) = heatmap(reshape(Array(x), Nx, Ny)', color=:viridis)
plotsol!(x, Nx=Nx, Ny=Ny; kwargs...) = heatmap!(reshape(Array(x), Nx, Ny)'; color=:viridis, kwargs...)

Nx = 151
	Ny = 100
	lx = 8pi
	ly = 2*2pi/sqrt(3)

function Laplacian2D(Nx, Ny, lx, ly, bc = :Neumann)
	hx = 2lx/Nx
	hy = 2ly/Ny
	D2x = CenteredDifference(2, 2, hx, Nx)
	D2y = CenteredDifference(2, 2, hy, Ny)
	if bc == :Dirichlet
		Qx = Dirichlet0BC(typeof(hx))
		Qy = Dirichlet0BC(typeof(hy))
	elseif bc == :Neumann
		Qx = Neumann0BC(hx)
		Qy = Neumann0BC(hy)
	elseif bc == :Periodic
		Qx = PeriodicBC(hx)
		Qy = PeriodicBC(hy)
	end
	# @show norm(D2x - Circulant(D2x[1,:]))
	A = kron(sparse(I, Ny, Ny), sparse(D2x * Qx)[1]) + kron(sparse(D2y * Qy)[1], sparse(I, Nx, Nx))
	return A, D2x
end

function F_sh(u, p)
	@unpack l, ν, L1 = p
	return -L1 * u .+ (l .* u .+ ν .* u.^2 .- u.^3)
end

function dF_sh(u, p)
	@unpack l, ν, L1 = p
	return -L1 .+ spdiagm(0 => l .+ 2 .* ν .* u .- 3 .* u.^2)
end

d2F_sh(u, p, dx1, dx2) = (2 .* p.ν .* dx2 .- 6 .* dx2 .* u) .* dx1
d3F_sh(u, p, dx1, dx2, dx3) = (-6 .* dx2 .* dx3) .* dx1
jet = (F_sh, dF_sh, d2F_sh, d3F_sh)

X = -lx .+ 2lx/(Nx) * collect(0:Nx-1)
Y = -ly .+ 2ly/(Ny) * collect(0:Ny-1)

sol0 = [(cos(x) .+ cos(x/2) * cos(sqrt(3) * y/2) ) for x in X, y in Y]
	sol0 .= sol0 .- minimum(vec(sol0))
	sol0 ./= maximum(vec(sol0))
	sol0 = sol0 .- 0.25
	sol0 .*= 1.7
	heatmap(sol0', color=:viridis)

Δ, D2x = Laplacian2D(Nx, Ny, lx, ly, :Neumann)
L1 = (I + Δ)^2
par = (l = -0.1, ν = 1.3, L1 = L1)

optnew = NewtonPar(verbose = true, tol = 1e-8, maxIter = 20)
# optnew = NewtonPar(verbose = true, tol = 1e-8, maxIter = 20, eigsolver = EigArpack(0.5, :LM))
	sol_hexa, hist, flag = @time newton(F_sh, dF_sh, vec(sol0), par, optnew)
	println("--> norm(sol) = ", norm(sol_hexa, Inf64))
	plotsol(sol_hexa)