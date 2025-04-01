export hold_segments!
export calculate_W

function calculate_W(p::EETCProblem, V::Real)
    W = p.train.ρ > 0. ? Roots.find_zero(v -> -ψ(p.train, V) + p.train.ρ*ψ(p.train, v), V) : V
end

function hold_segments!(p::EETCProblem, V::T) where {T<:Real}
    W = calculate_W(p, V)
    # Find points at which gradient and/or speed limit changes
    segment_borders = segmentize!(p.track)
    push!(segment_borders, p.track.length)

    ports = Port[]

    for (k, segstart) in enumerate(segment_borders[1:end-1])
    # Find at which points the train can hold speed V with u > 0
    # These segments fulfill U̅(V) - r(V) + g(x) > 0 && -r(V) + g(x) ≤ 0
        if p.train.U̅(V) - r(p.train, V) + g(p.track, segstart + 0.1) > 0 &&
            -r(p.train, V) + g(p.track, segstart) ≤ 0

            if speedlimit(p.track, segstart + 0.1) < V  # choose the lower of speed limit and V
                push!(ports, Port(segstart, segment_borders[k+1], HoldP_SL, speedlimit(p.track, segstart + 0.1)))
            else
                if !isempty(ports) && ports[end].speed == V && ports[end].mode == HoldP
                    ports[end] = Port(ports[end].start, segment_borders[k+1], HoldP, V)
                else
                    push!(ports, Port(segstart, segment_borders[k+1], HoldP, V))
                end
            end

        # Find at which points the train can hold speed W with u < 0
        # These segments fulfill -r(W) + g(x) ≥ 0
        elseif -r(p.train, W) + g(p.track, segstart + 0.1) ≥ 0

            if speedlimit(p.track, segstart + 0.1) < W
                push!(ports, Port(segstart, segment_borders[k+1], HoldR_SL, speedlimit(p.track, segstart + 0.1)))
            else
                push!(ports, Port(segstart, segment_borders[k+1], HoldR, W))
            end
        end
    end

    ports
end