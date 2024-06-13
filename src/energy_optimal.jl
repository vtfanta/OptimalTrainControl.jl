export hold_segments

function hold_segments(p::EETCProblem, V::T, W::S) where {T<:Real, S<:Real}
    # Find points at which gradient and/or speed limit changes
    segment_borders = segmentize!(p.track)
    push!(segment_borders, p.track.length)

    ports = Port[]

    for (k, segstart) in enumerate(segment_borders[1:end-1])
    # Find at which points the train can hold speed V with u > 0
    # These segments fulfill U̅(V) - r(V) + g(x) > 0 && -r(V) + g(x) ≤ 0
        if p.train.U̅(V) - r(p.train, V) + g(p.track, segstart + 0.1) > 0 &&
            -r(p.train, V) + g(p.track, segstart) ≤ 0

            if speedlimit(p.track, segstart + 0.1) < V
                push!(ports, Port(segstart, segment_borders[k+1], HoldP_SL, speedlimit(p.track, segstart + 0.1)))
            else
                push!(ports, Port(segstart, segment_borders[k+1], HoldP, V))
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