# This file is a part of SpacingStatistics.jl, licensed under the MIT License (MIT).

using SpacingStatistics
using Test


@testset "spacingstats" begin
    @test @inferred(SpacingStatistics.hello_world()) == 42
end
