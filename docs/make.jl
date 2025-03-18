using GBPlotsInteractive
using Documenter

DocMeta.setdocmeta!(GBPlotsInteractive, :DocTestSetup, :(using GBPlotsInteractive); recursive = true)

makedocs(;
    modules = [GBPlotsInteractive],
    authors = "jeffersonparil@gmail.com",
    sitename = "GBPlotsInteractive.jl",
    format = Documenter.HTML(;
        canonical = "https://GenomicBreeding.github.io/GBPlotsInteractive.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/GenomicBreeding/GBPlotsInteractive.jl", devbranch = "main")
