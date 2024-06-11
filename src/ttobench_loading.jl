using JSON

export load_ttobench_track

"""
    load_ttobench_track(filename; which_stops = (1, 2))

Load a track from a file in the TTOBench format.

# Arguments
- `filename::AbstractString`: path to the file containing the track.
- `which_stops::Tuple{Int, Int} = (1, 2)`: tuple of indices of the stops to consider (only journey between two stops is currently possible).
"""
function load_ttobench_track(filename::AbstractString; which_stops::Tuple{Int, Int} = (1, 2))
    if !isfile(filename)
        throw(ArgumentError("File $filename does not exist."))
    end

    top_dict = JSON.parsefile(filename)

    stops = top_dict["stops"]["values"]
    if Base.length(stops) < 2
        throw(ArgumentError("At least two stops are required."))
    end
    if !(stops[1] ≈ 0)
        throw(ArgumentError("First stop must be at position 0."))
    end

    @info "Loading track $filename from stop $(which_stops[1]) to stop $(which_stops[2])."

    stops = [stops[which_stops[1]], stops[which_stops[2]]]
    length = stops[which_stops[2]] - stops[which_stops[1]]
    length_unit = top_dict["stops"]["unit"]

    if length_unit == "km"
        length *= 1e3
        stops .*= 1e3
    elseif length_unit == "ft"
        length *= 0.3048
        stops .*= 0.3048
    elseif length_unit == "mi"
        length *= 1609.34
        stops .*= 1609.34
    end

    if "speed limits" ∉ keys(top_dict)
        throw(ArgumentError("Speed limits are not defined but are mandatory."))
    end

    # speed limits
    sl_position_unit = top_dict["speed limits"]["units"]["position"]
    sl_speed_unit = top_dict["speed limits"]["units"]["velocity"]

    sl_x = [sl[1] for sl in top_dict["speed limits"]["values"]]
    sl_v = [sl[2] for sl in top_dict["speed limits"]["values"]]

    if sl_position_unit == "km"
        sl_x .*= 1e3
    elseif sl_position_unit == "ft"
        sl_x .*= 0.3048
    elseif sl_position_unit == "mi"
        sl_x .*= 1609.34
    end

    if sl_speed_unit == "km/h"
        sl_v /= 3.6
    elseif sl_speed_unit == "ft/s"
        sl_v *= 0.3048
    elseif sl_speed_unit == "mph"
        sl_v *= 0.44704
    end

    track = Track(length = length, x_speedlimit = sl_x, speedlimit = sl_v)

    # altitude
    if "altitude" ∈ keys(top_dict)
        altitude = top_dict["altitude"]["value"]
        altitude_unit = top_dict["altitude"]["unit"]

        if altitude_unit == "km"
            altitude *= 1e3
        elseif altitude_unit == "ft"
            altitude *= 0.3048
        elseif altitude_unit == "mi"
            altitude *= 1609.34
        end

        track.altitude = altitude
    end

    # gradient
    if "gradients" ∈ keys(top_dict)
        g_position_unit = top_dict["gradients"]["units"]["position"]
        g_gradient_unit = top_dict["gradients"]["units"]["slope"]

        g_x = [g[1] for g in top_dict["gradients"]["values"]]
        g_g = [g[2] for g in top_dict["gradients"]["values"]]

        if g_position_unit == "km"
            g_x .*= 1e3
        elseif g_position_unit == "ft"
            g_x .*= 0.3048
        elseif g_position_unit == "mi"
            g_x .*= 1609.34
        end

        if g_gradient_unit == "permil"
            g_g /= 1000
        elseif g_gradient_unit == "percent"
            g_g /= 100
        end
            
        track = Track(length = length, x_speedlimit = sl_x, speedlimit = sl_v,
            x_gradient = g_x, gradient = g_g)
    end

    track
end