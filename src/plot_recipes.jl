using RecipesBase

export mode2color

# TODO add Vector{Mode} plotting recipe

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
    elseif m == HoldP_SL
        return :magenta
    elseif m == HoldR_SL
        return :purple
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

    if isnothing(track.x_gradient)
        X = [0.]
        @info "Plotting a flat track."
    else
        X = copy(track.x_gradient)
    end
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

@recipe function f(port::Port)
    color = mode2color(port.mode)

    @series begin
        seriestype := :path

        # Ignore in legends
        primary := false   

        linecolor --> color     
        [port.start, port.finish], [port.speed, port.speed]
    end
end

@recipe function f(ports::Vector{Port})
    X = []
    Y = []
    colors = []
    for port in ports
        push!(X, [port.start, port.finish])
        push!(Y, [port.speed, port.speed])
        push!(colors, mode2color(port.mode))
    end
    colors = hcat(colors...)

    @series begin
        seriestype := :path

        # Ignore in legends
        primary := false   

        linecolor --> colors
        linewidth --> 2
        X, Y
    end
end