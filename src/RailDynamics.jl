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

# export types
export Scenario, Vehicle, Tram, Train, Track, ControlLaw, Model, Resistance

abstract type Model end
abstract type Resistance end
abstract type Scenario end
abstract type Track end
abstract type Train <: Vehicle end
abstract type Tram <: Vehicle end
abstract type Vehicle end

include("Models.jl")
include("Tracks.jl")

end
