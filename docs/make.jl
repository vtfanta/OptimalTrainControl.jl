using OptimalTrainControl
using Documenter

DocMeta.setdocmeta!(OptimalTrainControl, :DocTestSetup, :(using OptimalTrainControl); recursive=true)

makedocs(;
    modules=[OptimalTrainControl],
    authors="Vít Fanta <fantavit@fel.cvut.cz>",
    repo="https://github.com/vtfanta/OptimalTrainControl.jl/blob/{commit}{path}#{line}",
    sitename="OptimalTrainControl.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vtfanta.github.io/OptimalTrainControl.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Problem Statement" => "problem.md",
        "Tutorial" => "tutorials.md",
        "Reference" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/vtfanta/OptimalTrainControl.jl",
    devbranch="main",
)
