using LinearAlgebra
import NonlinearSolve as NLS
import SimpleNonlinearSolve as SNLS
using OptimalTrainControl
using OrdinaryDiffEq
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

function find_ports_flat(track::Track, train::Train, V, v₀, vf)
    found_ports = Port[]

    # compute maximum possible holding speed
    max_hold_speed = Roots.find_zero(v -> train.U̅(v) - r(train, v), (0.0, 400*3.6))  # search for roots in the range of 0 to 400 km/h
    
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
                        HoldP_SL,
                        speedlim_speed
                    )
                else
                    speedlim_port = Port(
                        track.x_speedlimit[k],
                        track.x_speedlimit[k+1],
                        HoldP_SL,
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

V = 23.0
v₀ = 1.0
vf = 1.0

ports = find_ports_flat(track, train, V, v₀, vf)

## Connecting ports with optimal trajectory segments

# How to connect the ports with optimal trajectory segments?
# It is proven in Khmelnitsky, 2000 that two ports can be connected by at most 1 path (e.g. if port1 -> port2 -> port3 exists, then port1 -> port3 does not).
# It is also proven that there is only 1 combinded path connecting the start port to the finish port, which is the optimal trajectory.

# Build a graph of ports and find links (edges) between them.
using Graphs
graph = SimpleDiGraph(length(ports))

for port1_idx in eachindex(ports)
    for port2_idx in port1_idx+1:length(ports)
        port1, port2 = ports[port1_idx], ports[port2_idx]
        if !has_path(graph, port1_idx, port2_idx)   # possible new link

        end
    end
end

function calculate_W(train, V)
    if train.ρ > 0
        prob = Roots.ZeroProblem(v -> -ψ(train, V) + train.ρ * ψ(train, v), V)
        Roots.solve(prob)
    else
        V
    end
end

# First, we need to fix rhs_flat to properly handle parameters
function rhs_flat(states, p, x)
    t, v, η = states
    current_mode, train, constants, V, W = p.current_mode, p.train, p.costate_constants, p.V, p.W

    if current_mode == MaxP
        u = train.U̅(v)
        # costate η evaluation
        η_eq = (E(train, V, v) + last(constants)) / (u - r(train, v)) - η
    elseif current_mode == Coast
        u = zero(eltype(states))
        # costate η evaluation
        η_eq = (E(train, V, v) + last(constants)) / (u - r(train, v)) - η
    elseif current_mode == HoldP
        u = r(train, v)  # Fixed: was train.r(v)
        # costate η evaluation (calculate ζ and shift by ρ-1)
        η_eq = (train.ρ * E(train, W, v) + last(constants)) / (u - r(train, v)) + (train.ρ - 1) - η
    elseif current_mode == MaxB
        u = train.U̲(v)
        # costate η evaluation (calculate ζ and shift by ρ-1)
        η_eq = (train.ρ * E(train, W, v) + last(constants)) / (u - r(train, v)) + (train.ρ - 1) - η
    else
        error("Not implemented for mode: $(current_mode)")
    end

    return [
        1/v,
        (u - r(train, v)) / v,
        η_eq
    ]
end

flat_odefun = OrdinaryDiffEq.ODEFunction(rhs_flat, mass_matrix=LinearAlgebra.Diagonal([1, 1, 0]))

 @kwdef mutable struct EETCSimParams
    train
    costate_constants
    current_mode
    V    # optimal cruising speed
    W    # optimal braking speed
    track::OptimalTrainControl.Track
    ρ    # braking energy regeneration ratio ∈ [0, 1)
end

function define_callbacks(ρ)
    cb_lowspeed = ContinuousCallback(   # terminate at low speed to avoid singularity
        (states, params, x) -> states[2] - 1e-2,
        int -> OrdinaryDiffEq.SciMLBase.terminate!(int),
    )

    cb_maxp2coast = ContinuousCallback( # switch between MaxP and Coast when η crosses zero
        (states, params, x) -> states[3],
        int -> int.p.current_mode = MaxP,
        int -> int.p.current_mode = Coast
    )

    cb_coast2maxb = ContinuousCallback( # switch between MaxP and Coast when η crosses ρ - 1
        (states, params, x) -> states[3] - (ρ - 1.0),
        function affect_maxb2coast!(int)
            t, v, η = int.u
            x = int.t
            p = int.p

            p.current_mode = Coast

            # append constant since we're switching from ζ to η (from F to E type constant from Albrecht 2016)
            newE = E(p.train, p.V, v) - η * (0.0 - r(p.train, v) + OptimalTrainControl.g(p.track, x))
            push!(p.costate_constants, newE)
        end,
        function affect_coast2maxb!(int)
            t, v, η = int.u
            x = int.t
            p = int.p
            
            p.current_mode = MaxB  # change last(p.Es) to F from the Albrecht 2016 article

            # append to Es since switching from η to ζ (it's actually F)
            newF = (η - ρ + 1.0) * (p.train.U̲(v) - r(p.train, v) + OptimalTrainControl.g(p.track, x)) - ρ * E(p.train, p.W, v)
            push!(p.costate_constants, newF)
        end
    )

    cb_lowspeed, cb_maxp2coast, cb_coast2maxb
end

function get_init_E(mode::Mode, start_v, start_η, train, V, W)
    if mode == MaxP
        return start_η * (train.U̅(start_v) - r(train, start_v)) - E(train, V, start_v)
    elseif mode == Coast
        return start_η * (-r(train, start_v)) - E(train, V, start_v)
    elseif mode == MaxB
        return start_η * (train.U̲(start_v) - r(train, start_v)) - train.ρ * E(train, W, start_v)
    else
        error("Not implemented for mode: $(mode)")
    end
end

function simulate(xspan, train, init_states::Vector, start_mode::Mode, V, cb)
    W = calculate_W(train, V)
    init_constants = [get_init_E(start_mode, init_states[2], init_states[3], train, V, W)]  # initial value of the costate constant
    params = EETCSimParams(
        train=train,
        costate_constants=init_constants,
        current_mode=start_mode,
        V=V,
        W=W,
        track=track,
        ρ=train.ρ
    )
    cbs = define_callbacks(train.ρ)
    odeprob = ODEProblem{false}(rhs_flat, init_states, xspan, params)

    # Solve the ODE using a suitable solver
    odesol = OrdinaryDiffEq.solve(odeprob, Rodas5P(); dtmax=10.0, callback=CallbackSet(cbs..., cb))

    return odesol, params
end

##
track = Track(
    5e3;
    speedlimit=[40, 30, 13, 30., 40],
    x_speedlimit=[0, 1e3, 2e3, 3e3, 4e3]
)
V = 15.0
v₀ = 1.0
vf = 1.0
ports = find_ports_flat(track, train, V, v₀, vf)
PORT2_SPEED = V+30
PORT2_START = 0e3
PORT2_FINISH = 2e3
cb = ContinuousCallback(
    (states, params, x) -> states[2] - PORT2_SPEED,
    function (int) if PORT2_START ≤ int.t ≤ PORT2_FINISH
    terminate!(int) end
    end
)

# 5-element Vector{Port{Float64}}:
#  Port{Float64}(-Inf, 0.0, MaxP, 1.0)  P1
#  Port{Float64}(0.0, 2000.0, HoldP, 15.0)  P2
#  Port{Float64}(2000.0, 3000.0, HoldP_SL, 13.0)  P3
#  Port{Float64}(3000.0, 5000.0, HoldP, 15.0)  P4
#  Port{Float64}(5000.0, Inf, Coast, 1.0)  P5

#  P1 -> P2: MaxP to HoldP, 0 m -> 1307.69 m, found η(0) = 0.09206199978959455, within 1e-2 of the cruising speed, DONE
#  P2 -> P3: HoldP to HoldP_SL, not possible since transition to coasting at the point 1307.69 (previous port link) does not reach SL speed
#               in time, DONE
#  P1 -> P3: MaxP to HoldP_SL, 0 m -> 2000 m, need to find such η(0) that the trajectory reaches start of speedlimit port at the
#               speed limit speed. The η has discontinuity at the point of transition to HoldP_SL, such that it has value η(2000+) = 0.0
#               Found value: η(0) = 0.092005, within 1e-2 of the speed limit speed, DONE
#  P3 -> P4: HoldP_SL to HoldP with potential initial jump in η to MaxP to reach HoldP with η = 0 with the cruising speed
#               found: η(3000+) = 0.003082254 to reach V at 3480.273 m, DONE
#  P4 -> P5: HoldP to finish port, searching for location on P4 where to exit HoldP with η = 0 - eps() to coast and eventually MaxB
#               to reach end of track with final speed; not possible since transition to Coast at link from P3 does not reach final speed in time
#  P3 -> P5: HoldP_SL to finish port, not possible since even from start of P3 with η = 0, the trajectory does not reach final speed in time
#  P2 -> P5: HoldP to finish port; not possible since from the previous link from P1 to P2, the trajectory does not reach final speed in time
#  P1 -> P5: MaxP to finish port; finding η(0) such that the trajectory reaches final speed at the end of track; found η(0) = 0.0913696, within 1e-2 of the final speed
#                   found η(0) = 0.0913696, within 1e-2 of the final speed
#               NEED TO ALSO CHECK NO SPEEDLIMIT VIOLATIONS: there are none.
#                   if there were some, then there would be some η(0) such that the trajectory reaches the start of speedlimit at the appropriate
#                   speed limit speed, and then the trajectory would continue to the end of track, but there would be jump in η at the point of contact with the speed limit

odesol, params = simulate((.0, 5e3), train, [0.0, v₀, 0.09136959871495226], MaxP, V, cb); odesol

## root functions

function root_start_to_hold(η₀, port1, port2, train::Train, V, startmode::Mode)
    odesol, ret_params = simulate((port1.finish, port2.finish), train, [0.0, port1.speed, η₀], startmode, V, ContinuousCallback(
        (states, params, x) -> states[2] - port2.speed,
        function (int) if port2.start ≤ int.t ≤ port2.finish
            terminate!(int) end
        end
    ))
    if odesol.retcode == ReturnCode.Terminated
        if odesol[2,end] ≈ 1e-2     # terminated because speed too low
            return -Inf
        elseif abs(odesol[2,end] - port2.speed) > 1e-2   # terminated because speed did not reach V
            return sign(odesol[2,end] - port2.speed) * Inf
        else    # terminated because speed reached V
            return odesol[3,end] - 0.0   # η should be 0 at the point of reaching V, TODO how about speedhold
        end
    elseif odesol.retcode == ReturnCode.Success # reached end of port2
        return odesol[2,end] - port2.speed
    else
        error("ODE solver failed with retcode: $(odesol.retcode)")
    end
end

function root_start_to_finish(η₀, port1, port2, train::Train, V, startmode::Mode)
    odesol, ret_params = simulate((port1.finish, port2.start), train, [0.0, port1.speed, η₀], startmode, V, ContinuousCallback(
        (states, params, x) -> states[2] + 1e3, # never hit this callback
        int -> ()
    ))
    if odesol.retcode == ReturnCode.Terminated
        return odesol.t[end] - port2.start - port2.speed  # want to reach final speed at the end of track
    elseif odesol.retcode == ReturnCode.Success # reached finish
        return odesol[2,end] - port2.speed  # want to reach final speed at the end of track
    else
        error("ODE solver failed with retcode: $(odesol.retcode)")
    end
end

function root_start_to_speedlimit(η₀, port1, port2, train::Train, V, startmode::Mode)
    # find η(0) such that the trajectory reaches the start of speedlimit port at the speed limit speed
    odesol, _ = simulate((port1.finish, port2.start), train, [0.0, port1.speed, η₀], startmode, V, ContinuousCallback(
        (states, params, x) -> states[2] + 1e3, # never hit this callback
        int -> ()
    ))
    if odesol.retcode == ReturnCode.Terminated # exited early
        return odesol.t[end] - port2.start - port2.speed  # want to reach final speed at the end of track
    elseif odesol.retcode == ReturnCode.Success # reached start of speedlimit
        return odesol[2,end] - port2.speed  # want to reach final speed at the end of track
    else
        error("ODE solver failed with retcode: $(odesol.retcode)")
    end
end

function root_holdP_to_holdP_SL(x0, port1, port2, train::Train, V, startmode::Mode)
    # find x0 such that the trajectory reaches the start of speedlimit port at the speed limit speed
    # start with η=0 and coasting at x0
    odesol, _ = simulate((x0, port2.start), train, [0.0, port1.speed, 0.0-eps()], startmode, V, ContinuousCallback(
        (states, params, x) -> states[2] + 1e3, # never hit this callback
        int -> ()
    ))
    if odesol.retcode == ReturnCode.Terminated # exited early
        return odesol.t[end] - port2.start - port2.speed  # want to reach speedlimit speed at the start of port2
    elseif odesol.retcode == ReturnCode.Success # reached start of speedlimit
        return odesol[2,end] - port2.speed  # want to reach speedlimit speed at the start of port2
    else
        error("root_holdP_to_holdP_SL: ODE solver failed with retcode: $(odesol.retcode)")  
    end
end

function root_holdP_SL_to_subsequent_holdP(η_jump, port1, port2, train::Train, V, startmode::Mode)
    # find η_jump such that the trajectory reaches subsequent HoldP port at the cruising speed V with η=0
    odesol, _ = simulate((port1.finish, port2.finish), train, [0.0, port1.speed, 0.0+η_jump], startmode, V, ContinuousCallback(
        (states, params, x) -> states[2] - port2.speed,
        int -> ()
    ))
    if odesol.retcode == ReturnCode.Terminated # exited early
        return odesol[3,end] - 0.0   # want to reach η = 0 when V is reached
    elseif odesol.retcode == ReturnCode.Success # reached finish of subsequent HoldP
        return odesol[3,end]    # want to reach cruising speed at the start of port2
    else
        error("root_holdP_SL_to_subsequent_holdP: ODE solver failed with retcode: $(odesol.retcode)")  
    end
end

##
function link(port1, port2, train::Train, V; track=false)
    @assert port1.finish ≤ port2.start "Ports must be ordered: port1.finish ≤ port2.start"

    if port1.mode == HoldP_SL
        if port2.mode == HoldP && port1.finish ≈ port2.start    # connecting subsequent HoldP_SL to HoldP
            if port1.speed < port2.speed  # speedlimit speed less than cruising speed => need to accelerate to reach cruising speed
                # the above condition is always true, otherwise there is no need to connect the two ports.

                # Need to exit at port1.finish with v = port1.speed and I'm looking for jump in η such that port2 is reached with η = 0 and v = port2.speed
                η_jump_lower_bound = 0.0+eps()
                η_jump_upper_bound = 1e3

                sign_root_lower_bound = sign(root_holdP_SL_to_subsequent_holdP(η_jump_lower_bound, port1, port2, train, V, MaxP))
                sign_root_upper_bound = sign(root_holdP_SL_to_subsequent_holdP(η_jump_upper_bound, port1, port2, train, V, MaxP))

                if sign_root_lower_bound == sign_root_upper_bound
                    println("Rootfinding: Unable to find root for η_jump in root_holdP_SL_to_subsequent_holdP.")
                    return NaN
                end

                while root_holdP_SL_to_subsequent_holdP(η_jump_upper_bound, port1, port2, train, V, MaxP) > 0.0
                    η_jump_upper_bound /= 2.0
                end
                η_jump_upper_bound *= 2.0

                fn = NLS.NonlinearFunction((η,_) -> root_holdP_SL_to_subsequent_holdP(η[1], port1, port2, train, V, MaxP))
                root_prob = NLS.IntervalNonlinearProblem(fn, [η_jump_lower_bound, η_jump_upper_bound], [(η_jump_lower_bound+η_jump_upper_bound)/2], abstol=1e-6)
                return NLS.solve(root_prob).u
            end
        end

    elseif port1.mode == HoldP
        if port2.mode == HoldP_SL && port1.finish ≈ port2.start  # connecting subsequent HoldP to HoldP_SL
            if port1.speed > port2.speed  # cruising speed greater than upcoming speedlimit => need to slow down to reach start of speedlimit
                # the above condition is always true, otherwise there is no need to connect the two ports.
                # I can continue HoldP when there is no speedlimit requiring to slow down

                # need to find position such that speedlimit speed is achieved at the start of speedlimit port
                x0_lower_bound = port1.start + eps()
                x0_upper_bound = port1.finish - eps()

                sign_root_lower_bound = sign(root_holdP_to_holdP_SL(x0_lower_bound, port1, port2, train, V, Coast))
                sign_root_upper_bound = sign(root_holdP_to_holdP_SL(x0_upper_bound, port1, port2, train, V, Coast))
                if sign_root_lower_bound == sign_root_upper_bound
                    println("Rootfinding: Unable to find root for x0 in root_holdP_to_holdP_SL.")
                    return NaN
                end

                fn = NLS.NonlinearFunction((x,_) -> root_holdP_to_holdP_SL(x[1], port1, port2, train, V, Coast))
                root_prob = NLS.IntervalNonlinearProblem(fn, [x0_lower_bound, x0_upper_bound], [(x0_lower_bound+x0_upper_bound)/2], abstol=1e-2)
                return NLS.solve(root_prob).u

            else
                error("")
            end

        else
            error("Not implemented for connecting from HoldP to non-HoldP_SL port.")
        end

    elseif isinf(port1.start)
        if port2.mode == HoldP  # connecting from start port to HoldP port
            # find η₀ such that the trajectory reaches the start of HoldP port with η = 0 and speed = V
            if port1.mode == MaxP
                η₀_lower_bound = 0.0+eps()
                η₀_upper_bound = 0.1
                while root_start_to_hold(η₀_upper_bound, port1, port2, train, V, port1.mode) < 0.0
                    η₀_upper_bound *= 1.5
                    if η₀_upper_bound > 1e4
                        println("Rootfinding: Unable to find upper bound for η₀ in root_start_to_hold.")
                        return NaN
                    end
                end
                fn = NLS.NonlinearFunction((η,_) -> root_start_to_hold(η[1], port1, port2, train, V, port1.mode))
                root_prob = NLS.IntervalNonlinearProblem(fn, [η₀_lower_bound, η₀_upper_bound], [(η₀_lower_bound+η₀_upper_bound)/2], abstol=1e-6)
                return NLS.solve(root_prob).u
                # return Roots.find_zero(η₀ -> root_start_to_hold(η₀, port1, port2, train, V, port1.mode), 0.1, Order5(), atol=1e-6)
            elseif port1.mode == Coast
                η₀_upper_bound = 0.0-eps()
                η₀_lower_bound = -0.1
                while root_start_to_hold(η₀_lower_bound, port1, port2, train, V, port1.mode) > 0.0
                    η₀_lower_bound *= 1.5
                    if η₀_lower_bound < -1e4
                        println("Rootfinding: Unable to find lower bound for η₀ in root_start_to_hold.")
                        return NaN
                    end
                end
                fn = NLS.NonlinearFunction((η,_) -> root_start_to_hold(η[1], port1, port2, train, V, port1.mode))
                root_prob = NLS.IntervalNonlinearProblem(fn, [η₀_lower_bound, η₀_upper_bound], [(η₀_lower_bound+η₀_upper_bound)/2], abstol=1e-6)
                return NLS.solve(root_prob).u
                # return Roots.find_zero(η₀ -> root_start_to_hold(η₀, port1, port2, train, V, port1.mode), -0.1, Order5(), atol=1e-6)
            else
                error("Not implemented for starting mode: $(port1.mode)")
            end
        elseif isinf(port2.finish)  # connecting from start port to finish port
            # find η₀ such that the trajectory reaches the start of finish port with speed = vf
            if port1.mode == MaxP
                η₀_lower_bound = 0.0+eps()
                η₀_upper_bound = 0.1
                while root_start_to_finish(η₀_upper_bound, port1, port2, train, V, port1.mode) < 0.0
                    η₀_upper_bound *= 1.5
                    if η₀_upper_bound > 1e4
                        println("Rootfinding: Unable to find upper bound for η₀ in root_start_to_finish.")
                        return NaN
                    end
                end
                fn = NLS.NonlinearFunction((η,_) -> root_start_to_finish(η[1], port1, port2, train, V, port1.mode))
                root_prob = NLS.IntervalNonlinearProblem(fn, [η₀_lower_bound, η₀_upper_bound], [(η₀_lower_bound+η₀_upper_bound)/2], abstol=1e-6)
                return NLS.solve(root_prob).u
                # return Roots.find_zero(η₀ -> root_start_to_finish(η₀, port1, port2, train, V, port1.mode), (0.00000001, 10.0), Roots.Brent(), atol=1e-2)
            elseif port1.mode == Coast
                η₀_upper_bound = 0.0-eps()
                η₀_lower_bound = -0.1
                while root_start_to_finish(η₀_lower_bound, port1, port2, train, V, port1.mode) > 0.0
                    η₀_lower_bound *= 1.5
                    if η₀_lower_bound < -1e4
                        println("Rootfinding: Unable to find lower bound for η₀ in root_start_to_finish.")
                        return NaN
                    end
                end
                fn = NLS.NonlinearFunction((η,_) -> root_start_to_finish(η[1], port1, port2, train, V, port1.mode))
                root_prob = NLS.IntervalNonlinearProblem(fn, [η₀_lower_bound, η₀_upper_bound], [(η₀_lower_bound+η₀_upper_bound)/2], abstol=1e-6)
                return NLS.solve(root_prob).u
                # return Roots.find_zero(η₀ -> root_start_to_finish(η₀, port1, port2, train, V, port1.mode), (-10.0, -0.00000001), Roots.Brent(), atol=1e-2)
            else
                error("Not implemented for starting mode: $(port1.mode)")
            end
        elseif port2.mode == HoldP_SL  # connecting from start port to speedlimit port
            # find η₀ such that the trajectory reaches the start of speedlimit port at the speed limit speed
            if port1.mode == MaxP
                η₀_lower_bound = 0.0+eps()
                η₀_upper_bound = 0.1
                while root_start_to_speedlimit(η₀_upper_bound, port1, port2, train, V, port1.mode) < 0.0
                    η₀_upper_bound *= 1.5
                    if η₀_upper_bound > 1e4
                        println("Rootfinding: Unable to find upper bound for η₀ in root_start_to_speedlimit.")
                        return NaN
                    end
                end
                fn = NLS.NonlinearFunction((η,_) -> root_start_to_speedlimit(η[1], port1, port2, train, V, port1.mode))
                root_prob = NLS.IntervalNonlinearProblem(fn, [η₀_lower_bound, η₀_upper_bound], [(η₀_lower_bound+η₀_upper_bound)/2], abstol=1e-2)
                return NLS.solve(root_prob).u
                # return Roots.find_zero(η₀ -> root_start_to_speedlimit(η₀, port1, port2, train, V, port1.mode), 0.1, Roots.Order5(), atol=1e-2)
            elseif port1.mode == Coast
                η₀_upper_bound = 0.0-eps()
                η₀_lower_bound = -0.1
                while root_start_to_speedlimit(η₀_lower_bound, port1, port2, train, V, port1.mode) > 0.0
                    η₀_lower_bound *= 1.5
                    if η₀_lower_bound < -1e4
                        println("Rootfinding: Unable to find lower bound for η₀ in root_start_to_speedlimit.")
                        return NaN
                    end
                end
                fn = NLS.NonlinearFunction((η,_) -> root_start_to_speedlimit(η[1], port1, port2, train, V, port1.mode))
                root_prob = NLS.IntervalNonlinearProblem(fn, [η₀_lower_bound, η₀_upper_bound], [(η₀_lower_bound+η₀_upper_bound)/2], abstol=1e-2)
                return NLS.solve(root_prob).u
                # return Roots.find_zero(η₀ -> root_start_to_speedlimit(η₀, port1, port2, train, V, port1.mode), -0.1, Roots.Order5(), atol=1e-2)
            else
                error("Not implemented for starting mode: $(port1.mode)")
            end
        else
            error("Not implemented for connecting from start port to non-HoldP or non-finish port.")
        end
    else
        error("Not implemented for connecting from non-start port.")
    end
end