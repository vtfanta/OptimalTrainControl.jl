using OrdinaryDiffEq
using Roots

export solve

φ(train::Train, v) = v * r(train, v)
ψ(train::Train, v) = (train.r[2] + 2train.r[3]*v) * v^2
E(train::Train, V, v) = ψ(train, V)/v + OptimalTrainControl.r(train, v)

function _odefun_flat(s::A, p::EETCProblem, x::T) where {T<:Real, A<:AbstractArray{T,1}}
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
    ds2 = (u - r(p.train, v)) / v
    SA[ds1, ds2]
end

function __solve(p::EETCProblem{TV,S,U,Nothing,Nothing,VS}, V::A, V′::A) where {TV, S, U, VS, A<:AbstractFloat}
    # solve for given V and V'

    lowspeed_cond(s, x, int) = s[2] - 0.1
    lowspeed_aff(int) = terminate!(int)

    W = p.train.ρ > 0. ? Roots.find_zero(v -> -ψ(p.train, V) + p.train.ρ*ψ(p.train, v), V) : V

    ## MaxP phase
    targetspeed_cond(s, x, int) = s[2] - V′
    targetspeed_aff(int) = terminate!(int)
    targetspeed_cb = ContinuousCallback(targetspeed_cond, targetspeed_aff; affect_neg! = nothing)
    lowspeed_cb = ContinuousCallback(lowspeed_cond, lowspeed_aff)
    cbs = CallbackSet(targetspeed_cb, lowspeed_cb)

    maxprob = EETCProblem(p.T, p.train, p.track, MaxP)
    s0 = SA[0., p.initial_speed]
    xspan = (0., p.track.length)
    odeprob = ODEProblem(_odefun_flat, s0, xspan, maxprob)
    odesol = OrdinaryDiffEq.solve(odeprob, Tsit5(); callback = cbs)
    t1 = odesol[1,:]
    v1 = odesol[2,:]
    x1 = odesol.t
    η1 = (E.(p.train, V, v1) .- E(p.train, V, V′)) ./ (p.train.U̅.(v1) .- r.(p.train, v1))

    ## Coast phase
    L(v) = E(p.train, V, V) * (v - V) + φ(p.train, V)
    # v_coast2brake = 0.
    if V′ < V
        f(U) = V′ * (1. - (φ(p.train, V′) - p.train.ρ*φ(p.train, U))/
            (ψ(p.train, V) + φ(p.train, V′))) - U
        v_coast2brake = Roots.find_zero(f, V′)
    elseif p.train.ρ ≈ 0
        v_coast2brake = Roots.find_zero(v -> L(v), V)
    elseif p.train.ρ ≈ 1
        v_coast2brake = V
    else
        v_coast2brake = Roots.find_zero(v -> L(v) - p.train.ρ*φ(p.train,v), V)
    end
    
    targetspeed2_cond(s, x, int) = s[2] - v_coast2brake
    targetspeed2_aff(int) = terminate!(int)
    targetspeed2_cb = ContinuousCallback(targetspeed2_cond, targetspeed2_aff)
    cbs = CallbackSet(targetspeed2_cb, lowspeed_cb)

    coastprob = EETCProblem(p.T, p.train, p.track, Coast)
    s0 = SA[odesol[1,end], v1[end]]
    xspan = (x1[end], p.track.length)
    odeprob = ODEProblem(_odefun_flat, s0, xspan, coastprob)
    odesol = OrdinaryDiffEq.solve(odeprob, Tsit5(); callback = cbs)
    t2 = odesol[1,:]
    v2 = odesol[2,:]
    x2 = odesol.t
    η2 = (E.(p.train, V, v2) .- E(p.train, V, V′)) ./ (- r.(p.train, v2))

    ## MaxB phase
    targetspeed3_cond(s, x, int) = s[2] - 1.
    targetspeed3_aff(int) = terminate!(int)
    targetspeed3_cb = ContinuousCallback(targetspeed3_cond, targetspeed3_aff)
    cbs = CallbackSet(targetspeed3_cb, lowspeed_cb)

    brakeprob = EETCProblem(p.T, p.train, p.track, MaxB)
    s0 = SA[t2[end], v2[end]]
    xspan = (x2[end], p.track.length)
    odeprob = ODEProblem(_odefun_flat, s0, xspan, brakeprob)
    odesol = OrdinaryDiffEq.solve(odeprob, Tsit5(); callback = cbs)
    t3 = odesol[1,:]
    v3 = odesol[2,:]
    x3 = odesol.t
    η3 = p.train.ρ .- 1. .+
        (p.train.ρ .* E.(p.train, W, v3) .- p.train.ρ .* E.(p.train, W, v_coast2brake)) ./
        (p.train.U̲.(v3) .- r.(p.train, v3))

    return [(x1, t1, v1, η1), (x2, t2, v2, η2), (x3, t3, v3, η3)]
end

function _solve(p::EETCProblem{TV,S,U,Nothing,Nothing,VS}, V::A) where {TV, S, U, VS, A<:AbstractFloat}
    # solve for given V, insert HoldP if necessary and tune V' if necessary
    # and build OTCSolution

    V′ = V
    sols = __solve(p, V, V′)
    x1, t1, v1, η1 = sols[1]
    x2, t2, v2, η2 = sols[2]
    x3, t3, v3, η3 = sols[3]

    if x3[end] < p.track.length # x[end] < X
        Δx = p.track.length - x3[end]
        phases = [MaxP, HoldP, Coast, MaxB]
        x_phases = [x1[1], x1[end], x2[1] + Δx, x3[1] + Δx]

        Δt = Δx / V
        t2 .+= Δt
        t3 .+= Δt

        xs = [x1; x2 .+ Δx; x3 .+ Δx]
        ts = [t1; t2; t3]
        vs = [v1; v2; v3]
        ηs = [η1; η2; η3]

        s_tot = [[ts[k], vs[k]] for k in eachindex(xs)]

        odeprob = ODEProblem(_odefun_flat, SA[0., 1.], (0., p.track.length), p)
        odesol = DiffEqBase.build_solution(odeprob, Tsit5(), xs, s_tot, retcode = ReturnCode.Success)

        control = function (x)
            current_phase = phases[searchsortedlast(x_phases, x)]
            if current_phase == MaxP
                p.train.U̅(odesol(x)[2])
            elseif current_phase == HoldP
                r(p.train, odesol(x)[2])
            elseif current_phase == Coast
                zero(x)
            elseif current_phase == MaxB
                p.train.U̲(odesol(x)[2])
            end
        end

        sol_tot = OTCSolution(
            odesol,
            x_phases,
            phases,
            control,
            ηs)
        return sol_tot
    else    # x3[end] > X, need to adjust V'
        _fun = function(V′)
            sol = __solve(p, V, V′)
            x3 = sol[3][1]
            v3 = sol[3][3]
            if x3[end] < p.track.length
                return -Inf
            else x3[end] ≈ p.track.length
                return v3[end] - 1.
            end
        end
        
        opt_V′ = Roots.find_zero(_fun, [1.01, V])

        sols = __solve(p, V, opt_V′)
        x1, t1, v1, η1 = sols[1]
        x2, t2, v2, η2 = sols[2]
        x3, t3, v3, η3 = sols[3]

        xs = [x1; x2; x3]
        ts = [t1; t2; t3]
        vs = [v1; v2; v3]
        ηs = [η1; η2; η3]

        # @show _fun(opt_V′)

        s_tot = [[ts[k], vs[k]] for k in eachindex(xs)]

        odeprob = ODEProblem(_odefun_flat, SA[0., 1.], (0., p.track.length), p)
        odesol = DiffEqBase.build_solution(odeprob, Tsit5(), xs, s_tot, 
            retcode = ReturnCode.Success)

        phases = [MaxP, Coast, MaxB]
        x_phases = [x1[1], x2[1], x3[1]]

        control = function (x)
            current_phase = phases[searchsortedlast(x_phases, x)]
            if current_phase == MaxP
                p.train.U̅(odesol(x)[2])
            elseif current_phase == Coast
                zero(x)
            elseif current_phase == MaxB
                p.train.U̲(odesol(x)[2])
            end
        end

        sol_tot = OTCSolution(
            odesol,
            x_phases,
            phases,
            control,
            ηs)
        return sol_tot
    end

end

"""
    solve(problem::EETCProblem)

Compute `OTCSolution` of an energy-efficient train control `problem` on a flat track.

See also [`EETCProblem`](@ref), [`OTCSolution`](@ref).
"""
function solve(p::EETCProblem{TV,S,U,Nothing,Nothing,VS}; atol = 5.) where {TV,S,U,VS}
    # On a flat track, the mode sequence goes as MaxP -> (HoldP) -> Coast -> MaxB

    # Check feasibility of the problem
    TOsol = solve(TOTCProblem(p.train, p.track))
    # @show TOsol.odesol[1,end]
    if TOsol.odesol[1,end] ≥ p.T
        error("EETC problem infeasible due to the total time constraint ($(p.T) s).")
    end

    # first guess of the holding speed
    V = p.track.length / p.T 
    
    _f = function (V)
        sol = _solve(p, V)
        sol.odesol[1,end] - p.T
    end

    opt_V = Roots.find_zero(_f, V; atol)

    return _solve(p, opt_V)
end