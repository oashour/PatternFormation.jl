using PatternFormation
using Documenter

DocMeta.setdocmeta!(PatternFormation, :DocTestSetup, :(using PatternFormation); recursive=true)

makedocs(;
    modules=[PatternFormation],
    authors="Omar Ashour",
    repo="https://github.com/oashour/PatternFormation.jl/blob/{commit}{path}#{line}",
    sitename="PatternFormation.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://oashour.github.io/PatternFormation.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/oashour/PatternFormation.jl",
    devbranch="main",
)
