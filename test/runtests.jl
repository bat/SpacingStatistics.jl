# This file is a part of SpacingStatistics.jl, licensed under the MIT License (MIT).

import Test
import SpacingStatistics
import Documenter

Test.@testset "Package SpacingStatistics" begin
    include("test_spacingstats.jl")

    # doctests
    Documenter.DocMeta.setdocmeta!(
        SpacingStatistics,
        :DocTestSetup,
        :(using SpacingStatistics);
        recursive=true,
    )
    Documenter.doctest(SpacingStatistics)
end # testset
