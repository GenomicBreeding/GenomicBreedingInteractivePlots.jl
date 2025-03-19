module GBPlotsInteractive

using GBCore
using StatsBase, MultivariateStats, Distributions, LinearAlgebra, DataFrames, StatsBase
using UnicodePlots
using GLMakie, ColorSchemes

include("phenomes.jl")

export sparsity, removesparsestroworcol, loadphenomesdata, filterphenomesdata, addpc1pc2
export figurelayout, heatmapinteractive!, scatterplotinteractive!
export plotinteractive2d

end
