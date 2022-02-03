module PatternFormation

using StaticArrays
using LightLattices
using LoopVectorization
using FLoops
using Symbolics
using Distributions

export init_cond
export GS_Periodic!, GS_Neumann0!, NTF_Periodic!, Box!

include("InitCond.jl")
include("PDEDefs.jl")

end
