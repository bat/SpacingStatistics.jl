# This file is a part of SpacingStatistics.jl, licensed under the MIT License (MIT).

"""
    SpacingStatistics

Package containing spacing statistics.

RPS: Recursive Product of spacingstats
"""
module SpacingStatistics

using RelocatableFolders

using FileIO
using Interpolations
using Distributions
using ArgCheck
using CSV
using SpecialFunctions
using IntervalArithmetic
using Mmap
using SciPy
using Roots
using MappedArrays

export rps

const ASSETS = @path joinpath(@__DIR__, "../assets")

# include("spacingstats.jl")
include("recursive_product_spacings/recursive_product_spacings.jl")
include("sum_sorted_spacings/sum_sorted_spacings.jl")
include("product_complementary_spacings/product_complementary_spacings.jl")

end # module
