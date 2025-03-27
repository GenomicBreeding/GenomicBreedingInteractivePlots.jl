using Pkg
Pkg.add(url = "https://github.com/GenomicBreeding/GenomicBreedingCore.jl")
Pkg.develop(url = "https://github.com/GenomicBreeding/GenomicBreedingInteractivePlots.jl")
Pkg.add("Documenter")
Pkg.add("StatsBase")
Pkg.add("LinearAlgebra")
Pkg.add("Distributions")
using Documenter
using GenomicBreedingInteractivePlots

DocMeta.setdocmeta!(
    GenomicBreedingInteractivePlots,
    :DocTestSetup,
    :(using GenomicBreedingInteractivePlots);
    recursive = true,
)

makedocs(;
    modules = [GenomicBreedingInteractivePlots],
    authors = "jeffersonparil@gmail.com",
    sitename = "GenomicBreedingInteractivePlots.jl",
    format = Documenter.HTML(;
        canonical = "https://GenomicBreeding.github.io/GenomicBreedingInteractivePlots.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/GenomicBreeding/GenomicBreedingInteractivePlots.jl", devbranch = "main")
