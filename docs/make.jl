using Pkg
Pkg.add(url = "https://github.com/GenomicBreeding/GBCore.jl")
Pkg.develop(url = "https://github.com/GenomicBreeding/GBPlotsInteractive.jl")
Pkg.add("Documenter")
using Documenter
using GBCore
using GBPlotsInteractive
using StatsBase, MultivariateStats, Distributions, LinearAlgebra, DataFrames
using UnicodePlots
using GLMakie, ColorSchemes

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
