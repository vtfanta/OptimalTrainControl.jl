module Models

export AlbrechtModel

using Reexport
using Parameters

@reexport using RailDynamics

struct BasicScenario <: Scenario
    model::{<:Model}
    track::{<:Track}
    g<:Real
end

"""
Empirical formula originally calculated for freight cars. The resistance (in N/kg) is given by
    R = a + b * v + c * v^2,
where v is the vehicle speed.
"""
struct DavisResistance <: Resistance
    a::{<:Real}
    b::{<:Real}
    c::{<:Real}
end

struct AlbrechtModel <: Model
    resistance::{<:Resistance}
    maxcontrol::{<:Real}
    mincontrol::{<:Real}
    mass::{<:Real}
end

function odefun(s<:Scenario, u, p, t)
    m = s.model
    return function _odefun!(du, u, p, t)
        du[1] = inv(u[2])
        du[2] = (controllaw(m, u, p, t) - resistance(m.resistance, u) + inclinationforce(s, t)) * inv(u[2])
    end
end

function controllaw(m<:Model, u, p, t)
    # Pedal to the floor default
    m.maxcontrol
end

"""
Get the grade (in radians) of the track at the given position.
"""
function getgrade(track<:Track, t)
    # Flat default
    0
end

"""
Calculate the acting force component based on the current track grade.
The returned value is in N/kg, i.e. force per unit mass.
"""
function inclinationforce(s<:Scenario, t)
    - s.g * sin(getgrade(s.track, t))
end

"""
Calculates the Davis formula resistant force per unit mass.
    R = a + b * v + c * v^2,
where v is the vehicle speed and a, b and c are the resistance parameters.
"""
function resistance(r::DavisResistance, u)
    r.a + r.b * u[2] + r.c * u[2]^2
end