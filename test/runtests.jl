using Test
using Documenter
using GBCore
using GBPlotsInteractive
using StatsBase, MultivariateStats, Distributions, LinearAlgebra, DataFrames
using UnicodePlots
using GLMakie, ColorSchemes

Documenter.doctest(GBPlotsInteractive)

@testset "GBPlotsInteractive.jl" begin
    @test 1 == 1
end
