# This file is a part of SpacingStatistics.jl, licensed under the MIT License (MIT).

"""
    SpacingStatistics

Package containing spacing statistics.

RPS: Recursive Product of spacingstats
"""
module SpacingStatistics

using RelocatableFolders

export rps

const ASSETS = @path joinpath(@__DIR__, "../assets")

# include("spacingstats.jl")
include("rps/utils.jl")

end # module
