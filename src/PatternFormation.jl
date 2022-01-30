module PatternFormation

using StaticArrays
using LightLattices
using LoopVectorization
using FLoops
using Symbolics

export init_cond
export GS_Periodic!, GS_Neumann0!, GS_Neumann1!, GS_Neumann2!, GS_Periodic1!, GS_Periodic2!

include("InitCond.jl")
include("PDEDefs.jl")

end
