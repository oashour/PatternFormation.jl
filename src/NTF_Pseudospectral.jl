using ApproxFun
using LinearAlgebra
using Plots; gr()
using Distributions
using OrdinaryDiffEq

println("Fourier Stuff")
S = Fourier()
N = 512
x = points(S, N)
y = points(S, N)
D2x = Derivative(S,2)[1:n,1:n]
D2 = kron(D2x, D2x)
T = ApproxFun.plan_transform(S, N)
Ti = ApproxFun.plan_itransform(S, N)

println("Caches")
u = zeros(N,N)
v = similar(u)
u_n = similar(u)
v_n = similar(u)
û_n = similar(u)
v̂_n = similar(u)
û_1 = similar(u)
v̂_1 = similar(u)
u_d2 = zeros(N*N)
v_d2 = zeros(N*N)

println("Setting up problem parameters")
Mᵤ =  1e3
Mᵥ =  1.0
γᵤ =  1.0
γᵥ =  1.0
a = 0.5
c = 0.00125 

println("Initial Conditions")
r01 = 0.01 .+ rand(Uniform(-0.01/100, +0.01/100), N, N)
r02 = 0.25 .+ rand(Uniform(-0.25/100, +0.25/100), N, N)
r0 = cat(r01, r02, dims=3)

cache_ops = (D2,T,Ti, u, v, u_n, v_n, û_n, v̂_n, û_1, v̂_1, u_d2, v_d2)
p = (Mᵤ, Mᵥ, γᵤ, γᵥ, a, c, cache_ops)
function NTF_nl!(dr̂,r̂,p,t)
    Mᵤ, Mᵥ, γᵤ, γᵥ, a, c, cache_ops = p
    D2,T,Ti, u, v, u_n, v_n, û_n, v̂_n, û_1, v̂_1 = cache_ops
    û = @view r̂[:,:,1]
    v̂ = @view r̂[:,:,2]
    dû = @view dr̂[:,:,1]
    dv̂ = @view dr̂[:,:,2]
    mul!(u, Ti, û)
    mul!(v, Ti, v̂)

    @. u_n = 2u*(u-1)*(2u-1)
    mul!(û_n, T, u_n)
    mul!(u_d2, D2, vec(û))
    @. û_1 = (û_n-γᵤ*reshape(u_d2,N,N)) # Figure out how to do this so you don't have to reshape?
    mul!(û_1, D2, û_1)
    @. u_n = u*v
    mul!(u_n, T, u_n)
    dû = Mᵤ*û_1 - a*u_n + c
    
    @. v_n = 2v*(v-1)*(2v-1)
    mul!(v̂_n, T, v_n)
    mul!(v_d2, D2, vec(v̂))
    @. û_1 = (v̂_n-γᵥ*reshape(v_d2, N, N))
    mul!(v̂_1, D2, v̂_1)
    @. v_n = u*v #Make efficient, need only once
    mul!(v_n, T, v_n)
    dv̂ = Mᵥ*v̂_1 - a*v_n + c
end

using Symbolics
#println("Jacobian sparsity")
#dr0 = copy(r0)
#@time jac_sparsity = Symbolics.jacobian_sparsity((dr,r)->NTF_nl!(dr,r,p,0.0),dr0,r0)