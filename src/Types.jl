# export types
export Scenario, Vehicle, Tram, Train, Track, ControlLaw, Model, Resistance
export DavisResistance
export ControlModes
export ModelParams, NewModelParams, SolverParams
export TrainProblem
export Segment

abstract type Model end
abstract type Resistance end
abstract type Scenario end
abstract type Track end
abstract type Vehicle end
abstract type Train <: Vehicle end
abstract type Tram <: Vehicle end

ControlModes = Set([:MaxP, :HoldP, :HoldR, :Coast, :MaxB])

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

# To allow broadcasting
Base.broadcastable(r::DavisResistance) = Ref(r)

struct AlbrechtModel <: Model
    resistance::Resistance
    maxcontrol
    mincontrol
    mass::Real
    ρ::Real
end
AlbrechtModel(r, ma, mi, m) = AlbrechtModel(r, ma, mi, m, 0.0)

mutable struct ModelParams
    u
    r
    g
    ρ
    currentmode
end
Base.broadcastable(p::ModelParams) = Ref(p)

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
    b = find_zero((f, x -> derivative(f, x)), 10 * one(V), Roots.Newton())
    h(x) = b + m.ρ * ψ(m.resistance, x)
    W = find_zero((h, x -> derivative(h, x)), 10 * one(V), Roots.Newton())

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

mutable struct Segment
    start
    finish
    mode
    entry
    exit
end
Segment(s,f,m) = Segment(s,f,m,nothing,nothing)
function Base.show(io::IO, s::Segment) 
    show(io, s.start)
    print(io, " ~ ")
    printstyled(io, s.mode; color = s.mode == :HoldP ? :blue : :magenta)
    print(io, " ~ ")
    show(io, s.finish)
end
Base.broadcastable(s::Segment) = Ref(s)

struct NewModelParams
    umax
    umin
    resistance
    ρ
    track
    V
    vᵢ
    vf
end

mutable struct SolverParams
    modelparams::NewModelParams
    currentmode
end

@with_kw mutable struct TrainProblem
    umax = v -> 1 / max(5, v)
    umin = v -> -1 / max(5, v)
    resistance = DavisResistance(1e-2, 0, 1.5e-5)
    ρ = 0
    track
    T
    vᵢ = 1.0
    vf = 1.0
    states = nothing
    control = nothing
    switchingpoints = nothing
end