using Pkg
Pkg.activate(".")
try
    Pkg.add(url = "https://github.com/GenomicBreeding/GBCore.jl")
catch
    nothing
end
using GBCore
using StatsBase, MultivariateStats, Distributions, LinearAlgebra, DataFrames
using GLMakie, ColorSchemes
