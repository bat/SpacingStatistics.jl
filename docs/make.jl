# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [nonstrict] [fixdoctests]
#
# for local builds.

using Documenter
using SpacingStatistics

# Doctest setup
DocMeta.setdocmeta!(
    SpacingStatistics,
    :DocTestSetup,
    :(using SpacingStatistics);
    recursive=true,
)

makedocs(
    sitename = "SpacingStatistics",
    modules = [SpacingStatistics],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://bat.github.io/SpacingStatistics.jl/stable/"
    ),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    linkcheck = !("nonstrict" in ARGS),
    strict = !("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/bat/SpacingStatistics.jl.git",
    forcepush = true,
    push_preview = true,
)
