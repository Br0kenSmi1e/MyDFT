using MyDFT
using Documenter

DocMeta.setdocmeta!(MyDFT, :DocTestSetup, :(using MyDFT); recursive=true)

makedocs(;
    modules=[MyDFT],
    authors="Br0kenSmi1e",
    sitename="MyDFT.jl",
    format=Documenter.HTML(;
        canonical="https://Br0kenSmi1e.github.io/MyDFT.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Br0kenSmi1e/MyDFT.jl",
    devbranch="main",
)
