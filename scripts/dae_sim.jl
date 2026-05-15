# use Rodas4 method to solve DAE (EETC with two states and the η state as algebraic)

using LinearAlgebra
using NonlinearSolve
using OptimalTrainControl
using OrdinaryDiffEq
using Plots
using StaticArrays

module myU
export Max_u, Min_u
import Base.max

abstract type MaxControl end
abstract type MinControl end

(u::MaxControl)(v::Real) = error("Not implemented. Write method for your concrete type")
(u::MinControl)(v::Real) = error("Not implemented. Write method for your concrete type")

struct Max_u{T<:Real} <: MaxControl
    P::T
    v₀::T
end

(u::Max_u)(v::Real) = u.P / Base.max(u.v₀, v)

struct Min_u{T<:Real} <: MinControl
    Q::T
    v₀::T
end

(u::Min_u)(v::Real) = u.Q / Base.max(u.v₀, v)
end

module MyResistance
export DavisResistance, ψ, E

abstract type AbstractResistance end

(r::AbstractResistance)(v::Real) = error("Not implemented. Write method for your concrete type: $(typeof(r))")
ψ(r::AbstractResistance, v::Real) = error("Not implemented. Write method for your concrete type: $(typeof(r))")
E(r::AbstractResistance, v::Real) = error("Not implemented. Write method for your concrete type: $(typeof(r))")

struct DavisResistance{T<:Real} <: AbstractResistance
    a::T
    b::T
    c::T
end

(r::DavisResistance)(v::Real) = r.a + r.b * v + r.c * v^2

ψ(r::DavisResistance, v::Real) = v^2 * (r.b + 2r.c * v)
E(r::DavisResistance, V::Real, v::Real) = ψ(r, V) / v + r(v)
end

module MySim
    using ..myU
    using ..MyResistance
    using OptimalTrainControl
    using OrdinaryDiffEq
    using Roots
    using StaticArrays

    function _rhs(states, params, x)    # simulate DAE; states are time, speed, η (algebraic evolution)
        _, v, η = states
        r, umax, umin, V, W, current_mode, constants = params.r, params.u_max, params.u_min, params.V, params.W, params.current_mode, params.costate_constants
        track = params.track
        ρ = params.ρ
        if current_mode == MaxP
            u = umax(v)
            # costate η evaluation
            η_eq = (MyResistance.E(r, V, v) + last(constants)) / (u - r(v) + OptimalTrainControl.g(track, x)) - η
        elseif current_mode == Coast
            u = zero(typeof(v))
            # costate η evaluation
            η_eq = (MyResistance.E(r, V, v) + last(constants)) / (u - r(v) + OptimalTrainControl.g(track, x)) - η
        elseif current_mode == MaxB
            u = umin(v)
            # costate η evaluation (calculate ζ and shift by ρ-1)
            η_eq = (ρ * MyResistance.E(r, W, v) + last(constants)) / (u - r(v) + OptimalTrainControl.g(track, x)) - 1 + ρ - η
        end

        dt = 1 / v
        dv = (u - r(v) + OptimalTrainControl.g(track, x)) / v

        # SA[dt, dv, η_eq]  # this apparently does not work with DAE initializer
        @MVector [dt, dv, η_eq]

        # d_states[1] = 1/v
        # d_states[2] = (u - r(v)) / v
        # d_states[3] = (MyResistance.E(r, V, v) - MyResistance.E(r, V, V)) / (umax(v) - r(v)) - η

    end

    function simulate(states_0, params, pos_span)
        odeprob = ODEProblem(_rhs, states_0, pos_span, params)
        #TODO
    end

    mutable struct EETCSimParams{T<:Real, U⁺<:myU.MaxControl, U⁻<:myU.MinControl, R<:MyResistance.AbstractResistance}
        u_max::U⁺
        u_min::U⁻
        r::R
        costate_constants::Vector{T}
        current_mode::OptimalTrainControl.Mode
        V::T    # optimal cruising speed
        W::T    # optimal braking speed
        track::OptimalTrainControl.Track
        ρ::T    # braking energy regeneration ratio ∈ [0, 1)
    end

    function define_callbacks(ρ::T) where {T<:Real}
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
                newE = MyResistance.E(p.r, p.V, v) - η * (0.0 - p.r(v) + OptimalTrainControl.g(p.track, x))
                push!(p.costate_constants, newE)
            end,
            function affect_coast2maxb!(int)
                t, v, η = int.u
                x = int.t
                p = int.p
                
                p.current_mode = MaxB  # change last(p.Es) to F from the Albrecht 2016 article

                # append to Es since switching from η to ζ (it's actually F)
                newF = (η - p.ρ + 1.0) * (p.u_min(v) - p.r(v) + OptimalTrainControl.g(p.track, x)) - p.ρ * MyResistance.E(p.r, p.W, v)
                push!(p.costate_constants, newF)
            end
        )

        cb_lowspeed, cb_maxp2coast, cb_coast2maxb
    end
    
    function calculate_W(r::R, ρ::T, V::T) where {R<:MyResistance.AbstractResistance, T<:Real}
        if ρ > 0
            prob = Roots.ZeroProblem(v -> -MyResistance.ψ(r, V) + ρ * MyResistance.ψ(r, v), V)
            Roots.solve(prob)
        else
            V
        end
    end
end

##

function get_init_E(mode::OptimalTrainControl.Mode, start_x::T, start_v::T, start_η::T, simparams::MySim.EETCSimParams) where {T<:Real}
    if mode == MaxP
        start_η * (simparams.u_max(start_v) - simparams.r(start_v) + OptimalTrainControl.g(simparams.track, start_x)) - MyResistance.E(simparams.r, simparams.V, start_v)
    elseif mode == Coast
        start_η * (- simparams.r(start_v) + OptimalTrainControl.g(simparams.track, start_x)) - MyResistance.E(simparams.r, simparams.V, start_v)
    elseif mode == MaxB
        start_η * (simparams.u_min(start_v) - simparams.r(start_v) + OptimalTrainControl.g(simparams.track, start_x)) - MyResistance.E(simparams.r, simparams.V, start_v)
    end
end

using ..MySim
using Roots

M = SA[1. 0 0; 0 1 0; 0 0 0]

##

function linking(position::T) where T<:Real
    function rhs(states, params, x)
        _, v, η = states
        r, umax, umin, V, current_mode = params
        if current_mode == MaxP
            u = umax(v)
        elseif current_mode == Coast
            u = zero(typeof(v))
        elseif current_mode == MaxB
            u = umin(v)
        end

        dt = 1 / v
        dv = (u - r(v)) / v
        η_eq = (MyResistance.E(r, V, v) - MyResistance.E(r, V, V)) / (u - r(v)) - η

        SA[dt, dv, η_eq]

        # d_states[1] = 1/v
        # d_states[2] = (u - r(v)) / v
        # d_states[3] = (MyResistance.E(r, V, v) - MyResistance.E(r, V, V)) / (umax(v) - r(v)) - η

    end

    M = SA[1. 0 0; 0 1 0; 0 0 0]
    f = ODEFunction{false,SciMLBase.FullSpecialize}(rhs, mass_matrix=M)

    res = MyResistance.DavisResistance(1e-2, 0.0, 1.5e-5)
    u_max = myU.Max_u(1, 5)
    u_min = myU.Min_u(-1, 5)
    V = 25

    lowspeed_cb = ContinuousCallback(
        (states, params, x) -> states[2] - 1.0,
        int -> terminate!(int); abstol=0.2
    )

    # try to connect starting port with cruising part

    dist_span = reverse((0.0, position))
    # dist_span = reverse([0.0, 9e3])
    final_states = MArray{Tuple{3}}([0.0, V, 0.0])
    current_mode = MaxP
    params = (res, u_max, u_min, V, current_mode)
    prob = ODEProblem{false}(f, final_states, dist_span, params)

    sol = OrdinaryDiffEq.solve(prob, Rodas5(); callback=lowspeed_cb, save_everystep=false)

    if sol.retcode == ReturnCode.Success
        sol[2, end] - 1.0
    elseif sol.retcode == ReturnCode.Terminated
        -sol.t[end]  # want to hit x = 0
    else
        error("Unknown retcode behaviour.")
    end
end

##
function link(port1::Port{T}, port2::Port{T}, simparams::MySim.EETCSimParams{T, Uplus, Uminus, R}) where {T<:Real, Uplus, Uminus, R}
    if isinf(port1.start)   # port1 is starting point
        if isinf(port2.finish)  # try to directly connect start and finish points (rare)
            error("Not implemented.")
        end

        if port2.mode == HoldP
            # start from somewhere on port2 and simulate backwards to port1;
            # find root of linking function that returns negative numbers when terminated due to low speed
            # and positive when simulated to port1, but final(initial) speed is different

            f = ODEFunction{false,SciMLBase.FullSpecialize}(MySim._rhs, mass_matrix=M)
            init_η = 0.0
            init_states = SA[0.0, V, init_η]
            x_span = reverse((port1.finish, port2.finish))
            push!(simparams.costate_constants, get_init_E(simparams.current_mode, x_span[1], init_states[2], init_η, simparams))
            odeprob = ODEProblem{false}(f, init_states, x_span, simparams)

            cb_tuple = MySim.define_callbacks(simparams.ρ)

            cb_hit_singular = ContinuousCallback(
                (states, params, x) -> states[2] - simparams.V,   # hit v = V
                function affect!(int)
                    if 0.0 ≤ int.t ≤ 10e3   # in defined singular segment
                        SciMLBase.terminate!(int)
                    end
                end
            )

            function make_root_f(orig_prob::ODEProblem, cbs::CallbackSet)
                function root_f(new_x::T)   where {T<:Real}
                    newprob = remake(orig_prob; tspan = (new_x, orig_prob.tspan[2]))
                    sol::ODESolution = OrdinaryDiffEq.solve(newprob, OrdinaryDiffEq.Rodas5P(), callback=cbs, initializealg = SciMLBase.NoInit(),
                        save_everystep=false, save_start=false)
                    if sol.retcode == ReturnCode.Terminated
                        return newprob.tspan[2] - sol.t[end] - port1.speed    # - (init_speed), 1.0 in this case to make the function continuous
                    elseif sol.retcode == ReturnCode.Success
                        return sol[2,end] - port1.speed
                    else
                        error("Undefined behaviour for this return code.")
                    end
                end
            end

            my_f = make_root_f(odeprob, CallbackSet(cb_tuple...))

            zeroprob = ZeroProblem(my_f, (port2.finish + port2.start) / 2)

            x_root = Roots.solve(zeroprob, Steffensen(); atol=0.5)  # find such position that integrating backwards gives initial condition

            ret_prob = remake(odeprob; tspan=(x_root, odeprob.tspan[2]))
            odesol::ODESolution = OrdinaryDiffEq.solve(ret_prob, OrdinaryDiffEq.Rodas5P(), callback=CallbackSet(cb_tuple...), initializealg=SciMLBase.NoInit())
        end

    elseif isinf(port2.finish)  # port2 is finishing point
        # simulate forward from port1 to finish on port2.start (end of track) and matching port2.speed
            error("Not implemented.")
    else    # connect two interior ports (singular segments)
            error("Not implemented.")
    end

end
##

r = MyResistance.DavisResistance(1e-2, 0.0, 1.5e-5)
ρ = 0.5
V = 25.0
track = Track(30e3)
port_start = Port(-Inf, 0., MaxP, 1.0)
port_hold = Port(0., length(track), HoldP, 25.0)

simparams = MySim.EETCSimParams(
    myU.Max_u(1.0, 5.0),
    myU.Min_u(-1.0, 5.0),
    r,
    Float64[],
    MaxP,
    V,
    MySim.calculate_W(r, ρ, V),
    track,
    ρ
)

link(port_start, port_hold, simparams)