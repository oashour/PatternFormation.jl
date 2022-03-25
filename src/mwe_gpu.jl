using SparseArrays, CUDA, LinearAlgebra

CUDA.allowscalar(false)

N = 400
uA2 = cu(zeros(N,N))
A2u = cu(zeros(N,N))
A2 = cu(sprand(N, N, 0.01))
u = cu(rand(N, N))

mul!(A2u, A2, u)
mul!(uA2, u, A2')