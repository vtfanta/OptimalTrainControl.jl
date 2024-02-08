using RecipesBase

export mode2color

function mode2color(m::Mode)
    if m == MaxP
        return :green
    elseif m == MaxB
        return :red
    elseif m == Coast
        return :grey
    elseif m == HoldP
        return :blue
    elseif m == HoldR
        return :orange
    end
end

# Plot OTCSolution speed profile with color sections (modes)
@recipe function f(sol::OTCSolution)
    # xlabel --> "Distance (m)"
    # ylabel --> "Speed (m/s)"

    X = sol.odesol.t
    V = sol.odesol[2,:]

    linewidth --> 2

    color --> mode2color.(getindex(sol.phases,
        [searchsortedlast(sol.x_phases, x) for x in X]))
    
    X, V
end

# Plot track altitude profile
@recipe function f(track::Track)
    # xlabel --> "Distance (m)"
    # ylabel --> "Altitude (m)"

    X = copy(track.x_gradient)
    push!(X, track.length)

    Y = altitude.(track, X)
    minY = minimum(Y)
    maxY = maximum(Y)

    linecolor --> :black

    label --> false

    ylim --> (minY, maxY + 5(maxY - minY))

    alpha --> 0.5

    @series begin
        seriestype := :path

        # Ignore in legends
        primary := false   

        linecolor --> :black     
        fillcolor --> :lightgray
        fillrange := minY
        X, Y
    end
    X, Y
end
