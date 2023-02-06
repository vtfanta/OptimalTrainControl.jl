module RailDynamics


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
