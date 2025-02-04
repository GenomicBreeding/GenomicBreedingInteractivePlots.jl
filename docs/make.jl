using Pkg
Pkg.add(url="https://github.com/GenomicBreeding/GBCore.jl")
Pkg.develop(url="https://github.com/GenomicBreeding/GBPlots.jl")
Pkg.add("Documenter")
using GBPlots
using Documenter

DocMeta.setdocmeta!(GBPlots, :DocTestSetup, :(using GBPlots); recursive = true)

makedocs(;
    modules = [GBPlots],
    authors = "jeffersonparil@gmail.com",
    sitename = "GBPlots.jl",
    format = Documenter.HTML(;
        canonical = "https://GenomicBreeding.github.io/GBPlots.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/GenomicBreeding/GBPlots.jl", devbranch = "main")
