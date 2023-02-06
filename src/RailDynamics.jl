module RailDynamics

using Parameters

# export types
export Scenario, Vehicle, Tram, Train, Track, ControlLaw, Model, Resistance

abstract type Scenario end
abstract type Vehicle end
abstract type Tram <: Vehicle end
abstract type Train <: Vehicle end
abstract type Track end
abstract type Model end
abstract type Resistance end


"""
    hello()

Prints basic greeting.
"""
hello() = println("Hi, this is RaiajjjjjjjjjjlDynamics.jl! BROOOOOOOOO")

include("Models.jl")

end
