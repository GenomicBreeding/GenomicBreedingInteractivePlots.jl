using Test
using Documenter
using GenomicBreedingCore
using GenomicBreedingInteractivePlots
using StatsBase, MultivariateStats, Distributions, LinearAlgebra, DataFrames
using UnicodePlots
using GLMakie, ColorSchemes

Documenter.doctest(GenomicBreedingInteractivePlots)

@testset "GenomicBreedingInteractivePlots.jl" begin
    @test 1 == 1
end
