# module Models
using Reexport
using Parameters
using DifferentialEquations
using Unitful
using Roots
using Interpolations
using ForwardDiff: derivative

@reexport using RailDynamics

export BasicScenario, MinimalTimeScenario, OptimalScenario
export Resistance
export AlbrechtModel
export play, controllaw, calculatecontrol!

mutable struct OptimalScenario <: Scenario
    model::Model
    track::Track
    g::Real
    initialvalues
    finalvalues
    V
    W
    controllaw
end
function OptimalScenario(m, t, g, iv, fv, V)
    f(x) = x + ψ(m.resistance, V)
    b = find_zero((f, x -> derivative(f, x)), 10, Roots.Newton())
    h(x) = b + m.ρ * ψ(m.resistance, x)
    W = find_zero((h, x -> derivative(h, x)), 10, Roots.Newton())

    OptimalScenario(m, t, g, iv, fv, V, W, nothing)
end

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
    ρ::Real
end
AlbrechtModel(r, ma, mi, m) = AlbrechtModel(r, ma, mi, m, 0.0)

function play(s::BasicScenario, initialvalues)
    X = 100.0
    prob = ODEProblem(odefun(s), initialvalues, (0.0, X))
    solve(prob)
end

function play(s::MinimalTimeScenario)
    if isnothing(s.controllaw)
        error("The scenario does not have a control law. Try running `calculatecontrols!(scenario)`.")
    end

    prob = ODEProblem(odefun(s), s.initialvalues, (0.0, length(s.track)))
    sol = solve(prob)
end

function odefun(s::OptimalScenario)
    m = s.model
    return function (du, u, p, t)
        v = u[2]
        μ₂ = u[3]

        η = μ₂ * inv(v) - 1
        ζ = η + 1 - m.ρ

        μ₁ = -ψ(m.resistance, s.V)

        du[1] = inv(v)
        du[2] = (m.controllaw(t, v, η, ζ) - resistance(m.resistance, v) + 
            inclinationforce(s, t)) * inv(v)
        du[3] = μ₁ / v^2 + μ₂ * du[2] * inv(v) + 
            μ₂ * ψ(m.resistance, v) * inv(v)^3 - 
            (η > 0 ? η : 0) * derivative(m.maxcontrol, v) + 
            (ζ < 0 ? -ζ : 0) * derivative(m.mincontrol, v)
    end
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

function calculatecontrol!(s::OptimalScenario)
    # TODO
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
    brakeinterpolator = linear_interpolation(reverse(brakesol.t), reverse(brakesol.u))

    s.controllaw = 
        function (x) 
            vt = throttleinterpolator(x)
            vb = brakeinterpolator(x)
            if vt < vb - [0.1]
                s.model.maxcontrol(vt[1])
            else
                s.model.mincontrol(vb[1])
            end
        end
end

function controllaw(m::Model, u, p, t)
    # All gas, no brakes
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

function ψ(r::DavisResistance, u)
    if Base.length(u) > 1
        u[2]^2 * (r.b + 2r.c * u[2])
    else
        u[1]^2 * (r.b + 2r.c * u[1])
    end
end

# end # module