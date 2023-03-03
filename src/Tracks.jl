@reexport using RailDynamics

export FlatTrack, HillyTrack
export getgrade, inclinationforce, getgradientacceleration
export start, finish

# To allow broadcasting
Base.broadcastable(t::Track) = Ref(t)

struct FlatTrack <: Track 
    X::Real
end

struct HillyTrack <: Track 
    waypoints::DataFrame
    interpolator
end

HillyTrack(df::DataFrame) = HillyTrack(df, LinearInterpolator(df[:,:Distance], df[:,:Altitude], NoBoundaries()))

function HillyTrack(filepath::AbstractString) 
    dfcandidate = rename!(CSV.read(filepath, DataFrame), [:Distance, :Altitude])
    if length(names(dfcandidate)) > 2
        dfcandidate = rename!(CSV.read(filepath, DataFrame; transpose = true), 
            [:Distance, :Altitude])
    end
    HillyTrack(dfcandidate)
end

function length(t::FlatTrack)
    t.X
end

function length(t::HillyTrack)
    maximum(t.waypoints[:,:Distance])
end

function HillyTrack(X::AbstractVector{<:Real}, Y::AbstractVector{<:Real}) 
    df = DataFrame("Distance" => X, "Altitude" => Y) 
    HillyTrack(df, LinearInterpolator(X, Y, NoBoundaries()))
end

"""
Get the grade (in radians) of a flat track at the given position.
"""
function getgrade(track::Track, t)
    # Flat default track
    0
end

"""
    getgrade(track::HillyTrack, pos)

Get the track grade (in radians) of the `track` at the given position `pos`.
"""
function getgrade(track::HillyTrack, t)
    atan(gradient(track.interpolator, t))
end

"""
Calculate the acting force component based on the current track grade.
The returned value is in N/kg, i.e. force per unit mass.
"""
function inclinationforce(s::Scenario, t)
    @warn "This function is deprecated!"
    - s.g * sin(getgrade(s.track, t)) * inv(s.model.mass)
end

"""
    getgradientacceleration(::Scenario, position)

Return accelerating/decelerating component of the gravitational acceleration 
in meters per second squared [m s^-2] based on the gravitational acceleration value defined in the scenario. 
Positive means downhill, negative means uphill.
"""
function getgradientacceleration(s::Scenario, position)
    - s.g * sin(getgrade(s.track, position))
end

"""
    getgradientacceleration(::Track, position)

Return acceleration/decelerating component of the gravitational acceleration
`g = 9.81 m / s^2`. Positive means downhill, negative means uphill.
"""
function getgradientacceleration(t::Track, position)
    - 9.81 * sin(getgrade(t, position))
end

function gradient(int, t)
    (int(t + 0.1) - int(t - 0.1)) / 0.2
end

"""
    start(t::FlatTrack)

Return 0.0 m as the beginning of the flat track.

"""
function start(::FlatTrack)
    0.0
end

"""
    start(t::HillyTrack)

Return start of the track in metres.
"""
function start(t::HillyTrack)
    minimum(t.waypoints[!, :Distance])
end

"""
    finish(t::HillyTrack)

Return length of the track in metres.
"""
function finish(t::HillyTrack)
    maximum(t.waypoints[!, :Distance])
end

"""
    finish(t::FlatTrack)

Return length of the flat track in metres.
"""
function finish(t::FlatTrack)
    t.X
end

@recipe function f(track::FlatTrack)
    @info "You are plotting a flat track. You are going to get what you want."

    xlabel --> "Distance (m)"
    ylabel --> "Altitude (m)"

    linecolor --> :black

    label --> false

    @series begin
        seriestype := :path

        primary := false
        linecolor --> :black
        [start(track), finish(track)], [0.0, 0.0]
    end
end

@recipe function f(track::HillyTrack)
    xlabel --> "Distance (m)"
    ylabel --> "Altitude (m)"

    linecolor --> :black

    label --> false

    ylim --> (0.0, 5 * maximum(track.waypoints[!, :Altitude]))

    x = track.waypoints[!, :Distance]
    minalt = minimum(track.waypoints[!, :Altitude])
    if minalt ≤ 0.0
        y = track.waypoints[!, :Altitude] .+ minalt
    else
        y = track.waypoints[!, :Altitude]
    end

    @series begin
        seriestype := :path

        # Ignore in legends
        primary := false   
        linecolor --> :black     
        fillcolor --> :lightgray
        fillrange := 0 
        x, y
    end
    x, y
end