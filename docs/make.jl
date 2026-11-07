using StandingColumnWell
using Documenter

# Make `using StandingColumnWell` available to every doctest in docstrings and pages.
DocMeta.setdocmeta!(
    StandingColumnWell,
    :DocTestSetup,
    :(using StandingColumnWell);
    recursive = true,
)

makedocs(;
    modules = [StandingColumnWell],
    authors = "Gabriel-Dion <dion.gabriel100@gmail.com>",
    sitename = "StandingColumnWell.jl",
    format = Documenter.HTML(;
        canonical = "https://GHE-jl.github.io/StandingColumnWell.jl",
        edit_link = "main",
        assets = String[],
        mathengine = Documenter.KaTeX(),
        sidebar_sitename = false,
    ),
    pages = [
        "Home" => "index.md",
        "API reference" => "api.md",
    ],
    # Keep the build strict so broken cross-references or missing docstrings fail CI.
    checkdocs = :exports,
)

deploydocs(;
    repo = "github.com/GHE-jl/StandingColumnWell.jl",
    devbranch = "main",
)
