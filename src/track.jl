export r, g, isvalidposition, gradient, speedlimit, altitude, segmentize!, length

const GRAV_ACC = 9.81

"""
    length(track::Track)

Return length of `track` in metres.
"""
Base.length(t::Track) = t.length

"""
    r(train::Train, speed:Real)

Return resistance specific force of `train` at given `speed`.

The resistance specific force (per unit mass, i.e. acceleration) is positive and calculated as `train.r[1] + train.r[2] * v + train.r[3] * v^2`.
"""
@inline function r(train::Train, v::T) where {T<:Real}
    train.r[1] + train.r[2] * v + train.r[3] * v^2
end

"""
    g(track::Track, position::Real)

Return gravitational acceleration component at `position` regarding the gradient of `track`.
"""
@inline function g(track::Track, x::T) where {T<:Real}
    if !isvalidposition(track, x)
        throw(ArgumentError("Position $(x) out of bounds."))
    end

    -GRAV_ACC * sin(atan(gradient(track, x)))
end

"""
    isvalidposition(track::Track, position::Real)

Check if `position` is not out of bounds of the `track`.
"""
function isvalidposition(t::Track, position::T) where {T <: Real}
    if position < 0 || position > t.length
        return false
    end
    true
end

"""
    gradient(track::Track, position::Real)

Return gradient (in rise over run) at `position` on `track`.
For example, a gradient of a segment which rises 10 metres over 100 metres
of distance has gradient of 0.1.
"""
function gradient(t::Track, position::T) where {T <: Real}
    if !isvalidposition(t, position)
        throw(ArgumentError("Position $(position) out of bounds."))
    end
    isempty(t.x_gradient) && return(0.)
    
    return t.gradient[searchsortedlast(t.x_gradient, position)]
end

"""
    speedlimit(track::Track, position::Real)

Return speedlimit (in metres per second) at `position` on `track`.
"""
function speedlimit(t::Track, position::T) where {T <: Real}
    if !isvalidposition(t, position)
        throw(ArgumentError("Position $(position) out of bounds."))
    end
    isempty(t.speedlimit) && return Inf

    t.speedlimit[searchsortedlast(t.x_speedlimit, position)]
end

"""
    altitude(track::Track, position::Real)

Return altitude (in metres) at `position` on `track`.
"""
function altitude(t::Track, position::T) where {T <: Real}
    if !isvalidposition(t, position)
        throw(ArgumentError("Position $(position) out of bounds."))
    end
    if isempty(t.x_gradient)
        return t.altitude
    end

    integrator = t.altitude
    for k in 2:length(t.x_gradient)
        if t.x_gradient[k] > position
            midpoint = (t.x_gradient[k-1]+position)/2
            integrator += gradient(t, midpoint) * (position - t.x_gradient[k-1])
            break
        end
        midpoint = (t.x_gradient[k-1]+t.x_gradient[k])/2
        integrator += gradient(t, midpoint) * (t.x_gradient[k] - t.x_gradient[k-1])
    end
    if position > last(t.x_gradient)
        midpoint = (last(t.x_gradient)+position)/2
        integrator += gradient(t, midpoint) * (position - last(t.x_gradient))
    end
    return integrator
end

"""
    segmentize!(track::Track)

Modify `track.x_segments` such that its elements mark starts of
track parts on which both gradient and speed limit are constant.
"""
function segmentize!(t::Track) 
    if !isempty(t.x_segments)
        # already segmented
    elseif isempty(t.speedlimit) || isempty(t.speedlimit)
        t.x_segments = t.x_gradient
    else
        if !isempty(t.x_gradient)
            curr_x = t.x_gradient[1]
            @assert curr_x ≈ t.x_speedlimit[1] "First elements of `t.x_gradient` and t.x_speedlimit` do not align."

            t.x_segments = unique(sort([t.x_gradient; t.x_speedlimit]))
        else
            t.x_segments = t.x_speedlimit
        end
    end
end