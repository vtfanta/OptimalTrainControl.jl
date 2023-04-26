module RailDynamics

using BasicInterpolators
using CSV
using DataFrames
using DifferentialEquations
using ForwardDiff: derivative
using Parameters
using RecipesBase
using Reexport
using Roots
using Unitful

include("Types.jl")
include("Models.jl")
include("Tracks.jl")
include("Linking.jl")

end
