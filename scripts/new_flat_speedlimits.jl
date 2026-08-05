using OptimalTrainControl
using Roots

# How to solve EETC problem on a flat track with piecewise constant speedlimits?

# 1. Choose cruising speed V
# 2. Find possible ports.
#   If the cruising speed V is below or equal to the value of the maximum possible holding speed (root of U̅(v) - r(v) = 0), then
#   a cruising port spans the entire track (thanks to its flatness).
#   The other ports are the starting port, speedlimit holding ports and the finish port. Speedlimit ports are present on portions of the track where
#   the speed limit is below the cruising speed V.
#   The cruising ports with the active speed limit V̄
#   have the value of η corresponding to the active mode (cruising with positive u -> η=0, cruising with negative u -> η = ρ - 1.

function find_ports_flat(track::Track, train::Train, V::T, v₀::T, vf::T) where {T<:Real}
    found_ports = Port{T}[]

    # compute maximum possible holding speed
    max_hold_speed = Roots.find_zero(v -> train.U̅(v) - r(train, v), (0.0, 400*3.6))  # search for roots in the range of 0 to 400 km/h
    @debug "Maximum possible holding speed: $max_hold_speed"
    
    # starting port
    if v₀ > V   # coast from initial speed to cruising speed
        start_port = Port(-Inf, 0.0, Coast, v₀)
        push!(found_ports, start_port)
    elseif v₀ < V   # accelerate from initial speed to cruising speed
        start_port = Port(-Inf, 0.0, MaxP, v₀)
        push!(found_ports, start_port)
    elseif v₀ ≈ V   # no start_port, start cruising directly
    end

    # speedlimit ports
    if !isempty(track.speedlimit)
        # find speedlimit sections which are below the cruising speed V
        for (k, speedlim_speed) in enumerate(track.speedlimit)
            if speedlim_speed < V
                # add new speedlimit port

                if k > length(track.x_speedlimit) - 1
                    speedlim_port = Port(
                        track.x_speedlimit[k],
                        length(track),
                        HoldP,
                        speedlim_speed
                    )
                else
                    speedlim_port = Port(
                        track.x_speedlimit[k],
                        track.x_speedlimit[k+1],
                        HoldP,
                        speedlim_speed
                    )
                end

                push!(found_ports, speedlim_port)
            else # speedlim_speed >= V, no speedlimit port needed
                # add new cruising port
                if k > length(track.x_speedlimit) - 1
                    cruising_port = Port(
                        track.x_speedlimit[k],
                        length(track),
                        HoldP,
                        V
                    )
                else
                    cruising_port = Port(
                        track.x_speedlimit[k],
                        track.x_speedlimit[k+1],
                        HoldP,
                        V
                    )
                end
                push!(found_ports, cruising_port)
            end
        end
    else # no speedlimit ports
        if V < max_hold_speed   # able to hold cruising speed
            # add HoldP port spanning entire track
            hold_port = Port(0.0, length(track), HoldP, V)
            push!(found_ports, hold_port)
        end
    end

    # finish port 
    if vf > V   # accelerate to final speed from cruising speed 
        finish_port = Port(length(track), Inf, MaxP, vf)
        push!(found_ports, finish_port)
    elseif vf < V   # coast or coast/brake to final speed from cruising speed
        finish_port = Port(length(track), Inf, Coast, vf)
        push!(found_ports, finish_port)
    elseif vf ≈ V   # no finish_port, finish cruising directly
    end

    K = length(found_ports) - 1
    k = 1
    while k ≤ K
        port1, port2 = found_ports[k], found_ports[k+1]
        @show port1, port2
        if port1.mode == port2.mode && port1.speed == port2.speed
            # merge ports
            merged_port = Port(port1.start, port2.finish, port1.mode, port1.speed)
            found_ports[k] = merged_port
            deleteat!(found_ports, k+1)
            K -= 1
        else
            k += 1
        end
    end

    return found_ports
end

# test find_ports_flat
track = Track(
    5e3;
    speedlimit=[40, 30, 20, 30., 40],
    x_speedlimit=[0, 1e3, 2e3, 3e3, 4e3]
)

train = Train(v -> 1/v, v -> -1/v, (1e-2, 0., 1.5e-5), 0.3)

V = 49.0
v₀ = 1.0
vf = 1.0

find_ports_flat(track, train, V, v₀, vf)