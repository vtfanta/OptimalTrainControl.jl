using OrdinaryDiffEq
using Roots
using StaticArrays

lowspeed_cond(s, x, int) = s[2] - 0.5
lowspeed_aff(int) = terminate!(int)

function make_cond_coast2brake(train::Train, V::T) where {T<:Real}
    W = train.ρ > 0. ? Roots.find_zero(v -> -ψ(train, V) + train.ρ*ψ(train, v), V) : V
    return function cond_coast2brake(s, x, int)
        η = calculate_η(s, x, int, V, W)
        η - (int.p.train.ρ - 1.)
    end
end
aff_2coast!(int) = int.p.current_phase = Coast

function make_aff_2brake!(train::Train, V::T) where {T<:Real}
    W = train.ρ > 0. ? Roots.find_zero(v -> -ψ(train, V) + train.ρ*ψ(train, v), V) : V
    return function aff_2brake!(int)
        v, x = int.u[2], int.t
        if int.p.current_phase == Coast 
            int.p.current_phase = MaxB
            η = calculate_η(int.u, int.t, int, V, W)
            F = (η - int.p.train.ρ + 1.) * 
                (int.p.train.U̲(v) - r(int.p.train, v) + g(int.p.track, x))
                - int.p.train.ρ * E(int.p.train, W, v)
            push!(int.p.Es, F)
        else
            nothing
        end
    end
end

function calculate_η(s, x, int, V, W)
    v = s[2]
    
    if int.p.current_phase == MaxP
        val = (E(int.p.train, V, v) + last(int.p.Es)) / 
            (int.p.train.U̅(v) - OptimalTrainControl.r(int.p.train, v) + g(int.p.track, x))
        # @show val
    elseif int.p.current_phase == Coast
        val = (E(int.p.train, V, v) + last(int.p.Es)) / 
            (- OptimalTrainControl.r(int.p.train, v) + g(int.p.track, x))
    elseif int.p.current_phase == MaxB
        # have to calculate ζ, the Fs are also in the EETC.Es array
        val = (int.p.train.ρ*E(int.p.train, W, v) + last(int.p.Es)) /
            (int.p.train.U̲(v) - OptimalTrainControl.r(int.p.train, v) + g(int.p.train, x))
        # η = ζ + ρ - 1
        val += int.p.train.ρ - 1.
    end
    @show val
    return val
end

function _odefun(s::A, p::EETCProblem, x::T) where {T<:Real, A<:AbstractArray{T,1}}
    t, v = s

    if p.current_phase == MaxP
        u = p.train.U̅(v)
    elseif p.current_phase == Coast
        u = 0.
    elseif p.current_phase == MaxB
        u = p.train.U̲(v)
    end
    # @show u

    ds1 = 1/v
    ds2 = (u - r(p.train, v) + g(p.track, x)) / v
    SA[ds1, ds2]
end

function _link_start(prob::EETCProblem{TV,S,U,Nothing,Nothing,VS}, V::A, x0::A) where {TV, S, U, VS, A<:AbstractFloat}
    
    xspan = (x0, 0.)
    s0 = SA[0., V]
    
    p = EETCProblem(prob.T, prob.train, prob.track, MaxP)
    odeprob = ODEProblem(_odefun, s0, xspan, p)

    lowspeed_cb = ContinuousCallback(lowspeed_cond, lowspeed_aff)
    odesol = OrdinaryDiffEq.solve(odeprob, Tsit5(), 
        callback = lowspeed_cb)

    if odesol.retcode == ReturnCode.Terminated
        return -Inf
    elseif odesol.retcode == ReturnCode.Success
        return odesol[2,end] - prob.initial_speed
    end
end

function _link_finish(prob::EETCProblem{TV,S,U,Nothing,Nothing,VS}, V::A, x0::A) where {TV, S, U, VS, A<:AbstractFloat}
    
    xspan = (x0, prob.track.length)
    s0 = SA[0., V]
    
    p = EETCProblem(prob.T, prob.train, prob.track, Coast, V, [-E(prob.train, V, V)])
    odeprob = ODEProblem(_odefun, s0, xspan, p)

    lowspeed_cb = ContinuousCallback(lowspeed_cond, lowspeed_aff)
    coastbrake_cb = ContinuousCallback(make_cond_coast2brake(p.train, V), nothing,
        affect_neg! = make_aff_2brake!(p.train, V))

    cbs = CallbackSet(lowspeed_cb, coastbrake_cb)
    odesol = OrdinaryDiffEq.solve(odeprob, Tsit5(), 
        callback = cbs)

    if odesol.retcode == ReturnCode.Terminated
        return -Inf
    elseif odesol.retcode == ReturnCode.Success
        return odesol[2,end] - 1.0
    end
end

