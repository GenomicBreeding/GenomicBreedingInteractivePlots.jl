module GBPlotsInteractive

using GBCore
using StatsBase, MultivariateStats, Distributions, LinearAlgebra, DataFrames
using UnicodePlots
using GLMakie, ColorSchemes
# using PrecompileTools: @compile_workload

include("phenomes.jl")

export sparsity, removesparsestroworcol, loadphenomesdata, filterphenomesdata, addpc1pc2
export figurelayout, heatmapinteractive!, scatterplotinteractive!
export plotinteractive2d

end
