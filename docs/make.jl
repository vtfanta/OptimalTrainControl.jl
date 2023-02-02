using RailDynamics
using Documenter

DocMeta.setdocmeta!(RailDynamics, :DocTestSetup, :(using RailDynamics); recursive=true)

makedocs(;
    modules=[RailDynamics],
    authors="Vít Fanta <fantavit@fel.cvut.cz> and contributors",
    repo="https://github.com/vtfanta/RailDynamics.jl/blob/{commit}{path}#{line}",
    sitename="RailDynamics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vtfanta.github.io/RailDynamics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vtfanta/RailDynamics.jl",
    devbranch="main",
)
