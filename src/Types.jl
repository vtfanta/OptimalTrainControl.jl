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

ControlModes = Set([:MaxP, :HoldP, :HoldR, :Coast, :MaxB, :HoldPlim, :HoldRlim])

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

@with_kw mutable struct TrainProblem
    umax = v -> 3 / max(5, v)
    umin = v -> -3 / max(5, v)
    resistance = DavisResistance(1e-2, 0, 1.5e-5)
    ρ = 0
    track
    T
    vᵢ = 1.0
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