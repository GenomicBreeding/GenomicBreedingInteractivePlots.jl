using Pkg
Pkg.activate(".")
try
    Pkg.update()
catch
    nothing
end
using GenomicBreedingCore
using GenomicBreedingInteractivePlots
using StatsBase, MultivariateStats, Distributions, LinearAlgebra, DataFrames
using UnicodePlots
using GLMakie, ColorSchemes
