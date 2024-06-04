using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Hypatia

using FlagSOS, JLD2, JuMP, LinearAlgebra, ProgressMeter, CairoMakie

maxLvl = 6
solver = Hypatia.Optimizer

include("src/CornerCat456.jl")
include("src/TreeInducibility.jl")
include("src/BlockSizeTable.jl")
include("src/TreeProductTable.jl")
include("src/TreeProfiles.jl")
include("src/DrawProfiles.jl")
include("src/NonConvexity.jl")