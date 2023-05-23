# export types
export Scenario, Vehicle, Tram, Train, Track, ControlLaw, Model, Resistance
export DavisResistance
export ControlModes
export OldModelParams, ModelParams, SolverParams
export TrainProblem
export Segment

abstract type Model end
abstract type Resistance end
abstract type Scenario end
abstract type Track end
abstract type Vehicle end
abstract type Train <: Vehicle end
abstract type Tram <: Vehicle end

const ControlModes = Set([:MaxP, :HoldP, :HoldR, :Coast, :MaxB, :HoldPlim, :HoldRlim])

"""
    modecolor(mode ∈ ControlModes)

Return color corresponding to the `mode`.
"""
function modecolor(mode)
    if mode == :MaxP
        :green
    elseif mode == :HoldP
        :blue
    elseif mode == :HoldR
        :light_red
    elseif mode == :Coast
        :grey
    elseif mode == :MaxB
        :red
    elseif mode == :HoldPlim || mode == :HoldRlim
        :yellow
    end
end

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

mutable struct OldModelParams
    u
    r
    g
    ρ
    currentmode
end
Base.broadcastable(p::OldModelParams) = Ref(p)

"""
    Segment(start, finish, mode, holdspeed)

Define segment of a track between `start` and `finish` where `mode` is possible with
holding speed `holdspeed`.
"""
mutable struct Segment
    start
    finish
    mode
    holdspeed
end
Segment(s,f,m) = Segment(s,f,m,nothing)
function Base.show(io::IO, s::Segment) 
    show(io, s.start)
    print(io, " ~ ")
    printstyled(io, s.mode; color = modecolor(s.mode))
    print(io, " ~ ")
    show(io, s.finish)
end
Base.broadcastable(s::Segment) = Ref(s)

@with_kw mutable struct ModelParams
    umax = v -> 3 / max(5, v)
    umin = v -> -3 / max(5, v)
    resistance = DavisResistance(1e-2, 0, 1.5e-5)
    ρ = 0
    track
    V = 10
    vᵢ = 1.0
    vf = 1.0
    speedlimit = nothing
    speedlimitX = nothing
    speedlimitY = nothing
end

mutable struct SolverParams
    modelparams::ModelParams
    currentmode
end

"""
    TrainProblem(; track::Track, T, umax = v -> 3 / max(5, v), umin = v -> -3 / max(5, v),
        resistance = DavisResistance(1e-2, 0, 1.5e-5), ρ = 0, vᵢ = 1.0, vf = 1.0)

Define the optimal train control problem via specifying all of the relevant parameters as keyword arguments. All of the
arguments have default values except `track` and total journey time `T`.

# Arguments
- `track::Track`: Track specification
- `T`: Total time of journey in seconds
- `umax = v -> 3 / max(5, v)`: Maximum traction limit, output in (m/s^2)
- `umin = v -> -3 / max(5, v)`: Maximum braking traction limit, output in (m/s^2)
- `resistance = DavisResistance(1e-2, 0, 1.5e-5)`: Resistance function, output in (m/s^2)
- `ρ = 0`: Regeneration coefficient, 0 ≤ ρ < 1; tells how much braking energy can be reused
- `vᵢ = 1.0`: Initial speed (m/s)
- `vf = 1.0`: Final speed (m/s)
"""
@with_kw mutable struct TrainProblem
    "Maximum traction limit"
    umax = v -> 3 / max(5, v)
    "Maximum braking traction limit"
    umin = v -> -3 / max(5, v)
    "Resistance function"
    resistance = DavisResistance(1e-2, 0, 1.5e-5)
    "Regeneration coefficient (which portion of braking energy can be reused)"
    ρ = 0
    "Track gradient profile"
    track
    "Total journey time (seconds)"
    T
    "Initial speed (metres per second)"
    vᵢ = 1.0
    "Final speed (metres per second)"
    vf = 1.0
    speedlimit = nothing
    speedlimitX = nothing
    speedlimitY = nothing
    states = nothing
    control = nothing
    switchingpoints = nothing
end

function ModelParams(prob::TrainProblem, V)
    ModelParams(prob.umax, prob.umin, prob.resistance, prob.ρ, prob.track, 
        V, prob.vᵢ, prob.vf)
end