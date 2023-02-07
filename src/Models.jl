# module Models
using Reexport
using Parameters
using DifferentialEquations
using Unitful

@reexport using RailDynamics

export BasicScenario, DavisResistance, AlbrechtModel
export play, controllaw

mutable struct MinimalTimeScenario <: Scenario
    model::Model
    track::Track
    g::Real
    initialvalues
    finalvalues
    controllaw
end
MinimalTimeScenario(model, track, g, iv, fv) = MinimalTimeScenario(model, track, g, iv, fv, nothing)

struct BasicScenario <: Scenario
    model::Model
    track::Track
    g::Real
end
BasicScenario(m, t) = BasicScenario(m, t, 9.81u"m/s^2")

"""
Empirical formula originally calculated for freight cars. The resistance (in N/kg) is given by
    R = a + b * v + c * v^2,
where v is the vehicle speed.
"""
struct DavisResistance <: Resistance
    a::Real
    b::Real
    c::Real
end

struct AlbrechtModel <: Model
    resistance::Resistance
    maxcontrol::Real
    mincontrol::Real
    mass::Real
end

function play(s::BasicScenario, initialvalues)
    X = 100.0u"m"
    prob = ODEProblem(odefun(s), initialvalues, (0.0u"m", X))
    solve(prob)
end

function odefun(s::BasicScenario)
    m = s.model
    return function _odefun!(du, u, p, t)
        du[1] = inv(u[2])
        du[2] = (m.maxcontrol - resistance(m.resistance, u) + inclinationforce(s, t)) * inv(u[2])
    end
end

function odefun(s::MinimalTimeScenario)
    m = s.model
    function _odefun!(du, u, p, t)
        du[1] = inv(u[2])
        du[2] = (s.controllaw(u, p, t) - resistance(m.resistance, u) + 
            inclinationforce(s, t)) * inv(u[2])
    end
end

function calculatecontrol(s::MinimalTimeScenario)
    function _maxthrottle!(du, u, p, t)
        du[1] = inv(u[2])
        du[2] = (s.model.maxcontrol - resistance(m.resistance, u) + 
            inclinationforce(s, t)) * inv(u[2])
    end
    function _maxbrake!(du, u, p, t)
        du[1] = inv(u[2])
        du[2] = (s.model.mincontrol - resistance(m.resistance, u) + 
            inclinationforce(s, t)) * inv(u[2])
    end

    throttleprob = ODEProblem(_maxthrottle!, s.initialvalues, (0.0, length(s.track)))
    brakeprob = ODEProblem(_maxbrake!, s.finalvalues, (length(s.track), 0))

    
end

function controllaw(m::Model, u, p, t)
    # Pedal to the floor default
    m.maxcontrol
end

"""
Calculate the Davis formula resistant force per unit mass.
    R = a + b * v + c * v^2,
where v is the vehicle speed and a, b and c are the resistance parameters.
"""
function resistance(r::DavisResistance, u)
    r.a + r.b * u[2] + r.c * u[2]^2
end

# end # module