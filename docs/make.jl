using OptimalTrainControl
using Documenter

DocMeta.setdocmeta!(OptimalTrainControl, :DocTestSetup, :(using OptimalTrainControl); recursive=true)

makedocs(;
    modules=[OptimalTrainControl],
    authors="Vít Fanta <fantavit@fel.cvut.cz> and contributors",
    repo="https://github.com/vtfanta/OptimalTrainControl_v2.jl/blob/{commit}{path}#{line}",
    sitename="OptimalTrainControl.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Problem Statement" => "problem.md",
        "Tutorial" => "tutorial.md",
        "Reference" => "reference.md"
    ],
)

deploydocs(
    repo = "github.com/vtfanta/OptimalTrainControl_v2.jl.git",
)
