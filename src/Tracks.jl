using DataFrames
using CSV
using Interpolations
using Reexport

@reexport using RailDynamics

export FlatTrack, HillyTrack
export getgrade, inclinationforce

struct FlatTrack <: Track end

struct HillyTrack <: Track 
    waypoint::DataFrame
    interpolator
end

HillyTrack(df::DataFrame) = HillyTrack(df, linear_interpolation(df[:,:Distance], df[:,:Altitude];extrapolation_bc=Line()))

function HillyTrack(filepath::AbstractString) 
    dfcandidate = CSV.read(filepath, DataFrame)
    if length(names(dfcandidate)) > 2
        dfcandidate = CSV.read(filepath, DataFrame; transpose = true)
    end
    return rename(dfcandidate, [:Distance, :Altitude])
end

HillyTrack(X::AbstractVector{<:Real}, Y::AbstractVector{<:Real}) = HillyTrack(
    DataFrame("Distance" => X, "Altitude" => Y))

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
    atan(gradient(track.interpolator, t)[1])
end

"""
Calculate the acting force component based on the current track grade.
The returned value is in N/kg, i.e. force per unit mass.
"""
function inclinationforce(s::Scenario, t)
    - s.g * sin(getgrade(s.track, t)) * inv(s.model.mass)
end

