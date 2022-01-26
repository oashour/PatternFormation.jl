module PatternFormation

using StaticArrays
using LightLattices
using LoopVectorization
using FLoops

export init_cond
export GS_Periodic!, GS_Neumann0!

include("InitCond.jl")
include("PDEDefs.jl")

end
