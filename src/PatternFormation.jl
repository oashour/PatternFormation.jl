module PatternFormation

using StaticArrays
using LightLattices
using LoopVectorization
using FLoops
using Symbolics

export init_cond
export GS_Periodic!, GS_Neumann0!, GS_Test!

include("InitCond.jl")
include("PDEDefs.jl")

end
