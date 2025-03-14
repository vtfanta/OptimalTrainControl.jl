using CSV
using DataFrames
using LinearAlgebra
using OptimalTrainControl
using Plots

df = CSV.read("/home/vit/Downloads/TrackTest.csv", DataFrame)
track = Track(
    length = df.distances[end],
    x_gradient = df.distances[1:end-1],
    gradient = df.grades[1:end-1])