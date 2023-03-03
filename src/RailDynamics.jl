module RailDynamics

using Reexport
using Parameters
using DifferentialEquations
using Unitful
using Roots
using BasicInterpolators
using ForwardDiff: derivative
using DataFrames
using CSV
using BasicInterpolators
using Reexport
using RecipesBase

# export types
export Scenario, Vehicle, Tram, Train, Track, ControlLaw, Model, Resistance

abstract type Scenario end
abstract type Vehicle end
abstract type Tram <: Vehicle end
abstract type Train <: Vehicle end
abstract type Track end
abstract type Model end
abstract type Resistance end

include("Models.jl")
include("Tracks.jl")

end
