using OrdinaryDiffEq

export Mode, Train, Track, OTCSolution, TOTCProblem, EETCProblem, EETCSimParams
export MaxP, HoldP, HoldR, Coast, MaxB, HoldP_SL, HoldR_SL
export Port

"""
An Enum useful for specifying the five possible control modes:
- `MaxP`: Maximum power; usually takes place at the start.
- `HoldP`: Cruising; holding constant speed with positive input.
- `HoldR`: Regenerative braking; holding constant speed with negative input, can only take place on steep downhill segments.
- `Coast`: Coasting; neither applying power or braking.
- `MaxB`: Maximum braking; usually takes place at the end.
"""
@enum Mode begin
    MaxP = 0
    HoldP = 1
    HoldR = 2
    Coast = 4
    MaxB = 8
    HoldP_SL = 16
    HoldR_SL = 32
end    

"""
    train = Train(U̅, U̲, r, ρ = 0)

Defines a `train` to be used in `TOTCProblem` or `EETCProblem` construction

# Arguments
- `U̅::Function`: maximum traction specific force (per unit mass, i.e. acceleration) as a function of speed.
- `U̲::Function`: minimum traction specific force (braking, per unit mass, i.e. acceleration) as a function of speed.
- `r::Tuple`: triplet of coefficients defining resistance of the `train` as a quadratic function of the speed `resistance(v) = r[1] + r[2]v + r[3]v^2 `.
- `ρ::Real = 0`: regeneration coefficient ``ρ ∈ [0,1]`` specifying proportion of braking speed recovered (0 means no regeneration, 1 means all braking energy regenerated).

# Example
```julia
train = Train(
    U̅ = v ->  1/v,
    U̲ = v -> -1/v,
    r = (1e-2, 0., 1.5e-5),
    ρ = 0.3
)
```

See also [`TOTCProblem`](@ref), [`EETCProblem`](@ref).
"""
@kwdef struct Train{T<:Real,S<:Real,F1,F2}
    U̅::F1
    U̲::F2
    r::NTuple{3,T}
    ρ::S = zero(Float64)
end

# Train(U̅, U̲, r) = Train(U̅, U̲, r, 0)

# To allow broadcasting
Base.broadcastable(t::Train) = Ref(t)

"""
    track = Track(;length <keyword arguments>)

Defines a `track` to be used in `TOTCProblem` or `EETCProblem` construction



# Arguments
- `altitude::Real`: altitude of the start of the `track`.
- `x_gradient::Vector{Real}`: vector of positions at which the gradient changes.
- `gradient::Vector{Real}`: vector of grade values in rise/run, positive means uphill.
- `x_speedlimit::Vector{Real}`: vector of positions at which the speed limit changes.
- `speedlimit::Vector{Real}`: vector of speed limit values in ``\\mathrm{m/s}``.

# Examples
```julia
# Flat track
flat_track = Track(length = 1e3)
```

```julia
# Hilly track
hilly_track = Track(
    length = 3e3,
    altitude = 100.,
    x_gradient = [0.0, 1e3, 1.7e3],
    gradient = [2e-3, 0., 1e-3]
)
```

See also [`TOTCProblem`](@ref), [`EETCProblem`](@ref).
"""
# @kwdef mutable struct Track1{T<:Real, G<:Union{Nothing, Vector{T}},
#     S<:Union{Nothing, Vector{T}}}
#     length::T
#     altitude::T = 0.
#     x_gradient::G = nothing
#     gradient::G = nothing
#     x_speedlimit::S = nothing
#     speedlimit::S = nothing
#     x_segments::Union{Nothing,Vector{Any}} = nothing
# end

@kwdef mutable struct Track{T<:Real}
    length::T
    altitude::T
    x_gradient::Vector{T}
    gradient::Vector{T}
    x_speedlimit::Vector{T}
    speedlimit::Vector{T}
    x_segments::Vector{T}
end

function Track(length::T; 
    altitude::T = zero(T), 
    x_gradient::Vector{T} = T[],
    gradient::Vector{T} = T[],
    x_speedlimit::Vector{T} = T[],
    speedlimit::Vector{T} = T[],
    x_segments::Vector{T} = T[]) where {T<:Real}
Track(length, altitude, x_gradient, gradient, x_speedlimit, speedlimit, x_segments)
end

# Track(l, a, xg, g, xsl, sl) = Track(l, a, xg, g, xsl, sl, nothing)

# To allow broadcasting
Base.broadcastable(t::Track) = Ref(t)

struct ConcreteControlFunction{T<:Real,S<:Real,F1,F2,O<:OrdinaryDiffEq.ODESolution}
    x_mode_switch::Vector{T}
    modes_sequence::Vector{Mode}
    odesol::O
    train::Train{T,S,F1,F2}
end

function (u::ConcreteControlFunction)(x::T) where {T<:Real}
    current_phase = u.modes_sequence[searchsortedlast(u.x_mode_switch, x)]
    if current_phase == MaxP
        return (T)(u.train.U̅(odesol(x)[2]))
    elseif current_phase == Coast
        return zero(x)
    elseif current_phase == MaxB
        return (T)(u.train.U̲(odesol(x)[2]))
    end
end

function create_concrete_control_function(
    x_mode_switch::Vector{T},
    modes_sequence::Vector{Mode},
    odesol::OrdinaryDiffEq.ODESolution,
    train::Train{T,S,F1,F2}) where {T<:Real,S<:Real,F1,F2}

    ConcreteControlFunction(x_mode_switch, modes_sequence, odesol, train)
end

"""
    sol1::OTCSolution = solve(p1::TOTCProblem)
    sol2::OTCSolution = solve(p2::EETCProblem)

Is returned as a result of solving a `TOTCProblem` or a `EETCProblem`.

# Fields
- `odesol::SciMLBase.ODESolution`: solution of the differential equation coming from DifferentialEquations.jl package. The states are `[time, speed]` and are accessed with e.g. `sol1.u`. The distances along the track are accessed via e.g. `sol1.t`. The interpolation functionality should behave as expected although this is not guaranteed since the problem has multiple phases.
- `x_phases::Vector{Real}`: sequence of positions at which control modes (see [`Mode`](@ref)) changes.
- `phases::Vector{Mode}`: sequence of control modes (see [`Mode`](@ref)).
- `control::Function`: optimal control as a function of distance along the track.
- `η::Vector{Real}`: trajectory of the adjoint variable determining the current control mode (see [`Mode`](@ref)).

See also [`TOTCProblem`](@ref), [`EETCProblem`](@ref), [`Mode`](@ref).
"""
struct OTCSolution{T<:Real,S<:Real,F1,F2,O<:OrdinaryDiffEq.ODESolution}
    odesol::OrdinaryDiffEq.ODESolution
    x_phases::Vector{T}
    phases::Vector{Mode}
    control::ConcreteControlFunction{T,S,F1,F2,O}
    η::Vector{T}
end
OTCSolution(sol, x_phases, phases, control) = OTCSolution(sol, x_phases, phases, control, eltype(x_phases)[])

"""
    prob = TOTCProblem(;train::Train, track::Track, <keyword arguments>)

Formulate a time-optimal train control problem to be solved.

# Arguments
- `train::Train`: vehicle specification for the problem.
- `track::Track`: track specification for the problem
- `current_phase::Mode = MaxP`: (mainly for internal purposes) starting control mode.
- `initial_speed::Real = 1.`: starting speed; ``1 \\mathrm{m/s}``  is regarded as a stop.
"""
@kwdef mutable struct TOTCProblem{T,S,U,V<:AbstractFloat,F1,F2}
    train::Train{T,S,F1,F2}
    track::Track{U}
    current_phase::Mode = MaxP
    initial_speed::V = 1.
end

TOTCProblem(train, track) = TOTCProblem(train, track, MaxP, 1.)
TOTCProblem(train, track, mode) = TOTCProblem(train, track, mode, 1.)

"""
    prob = EETCProblem(; T, train, track, <keyword arguments>)

Formulate an energy-efficient train control problem.

# Arguments
- `T::Real`: total time of the trip.
- `train::Train`: train specification for the problem.
- `track::Track`: track specification for the problem.
- `current_phase::Mode = MaxP`: (mainly for internal purposes) starting control mode.
- `initial_speed::Real = 1.`: starting speed; ``1`` m/s is regarded as a stop.
- `Es::Vector{Real} = []`: (internal) vector of shifting constants used for calculation of the adjoint variable trajectory.
"""
# @kwdef mutable struct EETCProblem{TV,S,U,TG,TS,V<:AbstractFloat,W<:AbstractFloat}
#     T::V
#     train::Train{TV,S}
#     track::Track{U,TG,TS}
#     current_phase::Mode = MaxP
#     initial_speed::V = 1.
#     Es::Vector{W} = Float64[]
# end
@kwdef mutable struct EETCProblem{S<:Real,F1,F2}
    T::S
    train::Train{S,S,F1,F2}
    track::Track{S}
    initial_speed::S = one(S)
end

EETCProblem(T, train, track) = EETCProblem(T, train, track, one(typeof(T)))

@kwdef struct Port{T<:AbstractFloat}
    start::T
    finish::T
    mode::Mode
    speed::T
end

mutable struct EETCSimParams{S<:Real,F1,F2}
    eetcprob::EETCProblem{S,F1,F2}
    V::S
    W::S
    Es::Vector{S}
    current_phase::Mode
end