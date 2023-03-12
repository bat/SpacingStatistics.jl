# This file is a part of SpacingStatistics.jl, licensed under the MIT License (MIT).

using SpacingStatistics
using Test


@testset "spacingstats" begin
    events = [0.1, 0.4, 0.56, 0.3, 0.12]
    @test @inferred(SpacingStatistics.rps(events, Uniform())) == (0.901005261912951, 0.5912765565468424)
end
