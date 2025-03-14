using Pkg
Pkg.activate(".")
try
    Pkg.update()
catch
    nothing
end
using GBCore
using GBPlotsInteractive
using StatsBase, MultivariateStats, Distributions, LinearAlgebra, DataFrames
using UnicodePlots
using GLMakie, ColorSchemes
