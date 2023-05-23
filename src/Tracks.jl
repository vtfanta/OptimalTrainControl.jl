@reexport using OptimalTrainControl

export FlatTrack, HillyTrack
export getgrade, inclinationforce, getgradientacceleration
export start, finish, setspeedlimits!, subtrack

function setspeedlimits!(prob::TrainProblem, Xs, limits)
    @unpack track = prob
    @assert Base.length(limits) == Base.length(Xs) + 1 "Wrong dimensions. Xs define points of change of the speedlimit."
    
    @assert issorted(Xs)
    # @assert Xs[1] > start(track) && Xs[end] < finish(track)

    Xs[1] ≈ start(track) ? nothing : pushfirst!(Xs, start(track))
    Xs[end] ≈ finish(track) ? nothing : push!(Xs, finish(track))

    function speedlim(x)
        idx = BasicInterpolators.findcell(x, Xs, Base.length(Xs))
        return limits[idx]
    end
    prob.speedlimit = speedlim
    prob.speedlimitX = Xs
    prob.speedlimitY = limits
end

function setspeedlimits!(params::ModelParams, Xs, limits)
    @unpack track = params
    @assert Base.length(limits) == Base.length(Xs) + 1 "Wrong dimensions. Xs define points of change of the speedlimit."
    
    @assert issorted(Xs)
    # @assert Xs[1] > start(track) && Xs[end] < finish(track)

    Xs[1] ≈ start(track) ? nothing : pushfirst!(Xs, start(track))
    Xs[end] ≈ finish(track) ? nothing : push!(Xs, finish(track))

    function speedlim(x)
        idx = BasicInterpolators.findcell(x, Xs, Base.length(Xs))
        return limits[idx]
    end
    params.speedlimit = speedlim
    params.speedlimitX = Xs
    params.speedlimitY = limits
end

# To allow broadcasting
Base.broadcastable(t::Track) = Ref(t)

"""
    FlatTrack(X::Real)

Create instance of a flat track of length `X`.
"""
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
    if Base.length(names(dfcandidate)) > 2
        dfcandidate = rename!(CSV.read(filepath, DataFrame; transpose = true), 
            [:Distance, :Altitude])
    end
    HillyTrack(dfcandidate)
end

function length(t::FlatTrack)
    t.X
end

function length(t::HillyTrack)
    maximum(t.waypoints[:,:Distance]) - minimum(t.waypoints[:,:Distance])
end

"""
    HillyTrack(Xs, Ys)

Create instance of a track. The track is defined by a set of waypoints
along distance `Xs` and heights `Ys` at the appropriate x-points.
"""
function HillyTrack(X::AbstractVector, Y::AbstractVector) 
    f(xprev, xnow, yprev, γ) = tan(asin(γ / -9.81))*(xnow - xprev) + yprev
    function getys(γs, X)
        @assert Base.length(X) == Base.length(γs) + 1 "Lengths don't match!"
        ret = []
        for idx in eachindex(X)
            if idx == 1
                push!(ret, 0.0)
            else
                push!(ret, f(X[idx-1],X[idx],ret[idx-1],γs[idx-1]))
            end
        end
        return ret
    end
    if Base.length(X) == Base.length(Y)
        df = DataFrame("Distance" => X, "Altitude" => Y) 
        HillyTrack(df, LinearInterpolator(X, Y, NoBoundaries()))
    elseif Base.length(X) - 1 == Base.length(Y) # Y vector consists of gradients
        newY = getys(Y, X)
        HillyTrack(X, newY)
    else
        error("Wrong dimensions.")
    end
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

Return 0 as the beginning of the flat track `t`.

"""
function start(::FlatTrack)
    0.0
end

"""
    start(t::HillyTrack)

Return start of the track `t` in metres.
"""
function start(t::HillyTrack)
    minimum(t.waypoints[!, :Distance])
end

"""
    finish(t::HillyTrack)

Return length of the track `t` in metres.
"""
function finish(t::HillyTrack)
    maximum(t.waypoints[!, :Distance])
end

"""
    finish(t::FlatTrack)

Return length of the flat track `t` in metres.
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

    ylim --> (0.0, 5 * maximum(abs.(track.waypoints[!, :Altitude])))

    x = track.waypoints[!, :Distance]
    minalt = minimum(track.waypoints[!, :Altitude])
    if minalt ≤ 0.0
        y = track.waypoints[!, :Altitude] .- minalt
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

"""
    my_track(x)

Return height value at the distance `x`.
"""
function (t::HillyTrack)(x)
    t.interpolator(x)
end

"""
    my_flattrack(x)

Return 0 height value at the distance `x`.
"""
function (t::FlatTrack)(x)
    0
end

"""
    subtrack(full_track::HillyTrack, from, to)

Create subtrack of type `HillyTrack` of `full_track` spanning from `from` to `to`.
"""
function subtrack(t::HillyTrack, from, to)
    
    if from in t.waypoints[:,"Distance"]
        idx = findfirst(t.waypoints[:,"Distance"] .== from)
        startpoint = (t.waypoints[idx,"Distance"], t.waypoints[idx,"Altitude"])
    else
        startpoint = (from, t(from))
    end

    if to in t.waypoints[:,"Distance"]
        idx = findfirst(t.waypoints[:,"Distance"] .== to)
        endpoint = (t.waypoints[idx,"Distance"], t.waypoints[idx,"Altitude"])
    else
        endpoint = (to, t(to))
    end

    newwaypointsX = t.waypoints[from .< t.waypoints[:,"Distance"] .< to,"Distance"]
    newwaypointsY = t.waypoints[from .< t.waypoints[:,"Distance"] .< to,"Altitude"]

    HillyTrack([startpoint[1]; newwaypointsX; endpoint[1]],
        [startpoint[2]; newwaypointsY; endpoint[2]])
end

"""
    subtrack(full_track::HillyTrack, from)

Create subtrack of type `HillyTrack` of `full_track` from `from` until `finish(full_track)`.
"""
function subtrack(t::HillyTrack, from)
    subtrack(t, from, finish(t))
end

""" 
    subtrack(full_flattrack::FlatTrack, from)

Create subtrack of type `FlatTrack` of `full_flattrack` from `from` until `finish(full_flattrack)`.
"""
function subtrack(t::FlatTrack, from)
    FlatTrack(t.X - from)
end