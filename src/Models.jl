# module Models
using Reexport
using Parameters
using DifferentialEquations
using Unitful
using Interpolations

@reexport using RailDynamics

export BasicScenario, DavisResistance, AlbrechtModel, MinimalTimeScenario
export play, controllaw, calculatecontrol!

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
    maxcontrol
    mincontrol
    mass::Real
end

function play(s::BasicScenario, initialvalues)
    X = 100.0
    prob = ODEProblem(odefun(s), initialvalues, (0.0, X))
    solve(prob)
end

function play(s::MinimalTimeScenario)
    prob = ODEProblem(odefun(s), s.initialvalues, (0.0, length(s.track)))
    sol = solve(prob)
end

function odefun(s::MinimalTimeScenario)
    m = s.model
    return function (du, u, p, t)
        du[1] = (s.controllaw(t) - resistance(m.resistance, u[1]) + inclinationforce(s, t)) * inv(u[1])
    end
end

function odefun(s::BasicScenario)
    m = s.model
    return function (du, u, p, t)
        du[1] = inv(u[2])
        du[2] = (m.maxcontrol(u[2]) - resistance(m.resistance, u) + inclinationforce(s, t)) * inv(u[2])
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

function calculatecontrol!(s::MinimalTimeScenario)
    function _maxthrottle!(du, u, p, t)
        du[1] = (s.model.maxcontrol(u[1]) - resistance(s.model.resistance, u) + 
            inclinationforce(s, t)) * inv(u[1])
    end
    function _maxbrake!(du, u, p, t)
        du[1] = (s.model.mincontrol(u[1]) - resistance(s.model.resistance, u) + 
            inclinationforce(s, t)) * inv(u[1])
    end

    throttleprob = ODEProblem(_maxthrottle!, s.initialvalues, (0.0, length(s.track)))
    brakeprob = ODEProblem(_maxbrake!, s.finalvalues, (length(s.track), 0))

    throttlesol = solve(throttleprob)
    brakesol = solve(brakeprob)

    throttleinterpolator = linear_interpolation(throttlesol.t, throttlesol.u)
    brakeinterpolator = linear_interpolation(reverse(cat(brakesol.t, [0.0], dims = 1)), reverse(cat(brakesol.u, [last(brakesol.u)], dims = 1)))

    s.controllaw = 
        function (v) 
            vt = throttleinterpolator(v)
            vb = brakeinterpolator(v)
            if vt ≤ vb
                s.model.maxcontrol(vt[1])
            else
                s.model.mincontrol(vb[1])
            end
        end
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
    if Base.length(u) > 1
        r.a + r.b * u[2] + r.c * u[2]^2
    else
        r.a + r.b * u[1] + r.c * u[1]^2
    end
end

# end # module