using DifferentialEquations
import NumericalIntegration: integrate
using Parameters
using Plots
using OptimalTrainControl
using Roots
import ForwardDiff: derivative

function odefun!(du, u, p::SolverParams, x)
    t, v, η = u

    V = p.modelparams.V
    ρ = p.modelparams.ρ

    du[1] = inv(v)
    du[2] = (control(v, x, p) - resistance(p.modelparams.resistance, v) + 
        getgradientacceleration(p.modelparams.track, x)) * inv(v)

    if p.currentmode != :MaxB
        du[3] = (ψ(p.modelparams.resistance, v) - v^2 * derivative(v -> control(v, x, p), v)) * η / v^3 + (ψ(p.modelparams.resistance,v) - ψ(p.modelparams.resistance,V))/v^3
    else
        du[3] = (ψ(p.modelparams.resistance, v) - v^2 * derivative(v -> control(v, x, p), v)) * η / v^3 + (ψ(p.modelparams.resistance,v) - ψ(p.modelparams.resistance,V))/v^3 - (1 - ρ) * derivative(v -> control(v, x, p), v) / v
    end
end

function control(v, x, params::SolverParams)
    mode = params.currentmode
    if mode == :MaxP
        params.modelparams.umax(v)
    elseif mode == :HoldP || mode == :HoldR
        resistance(params.modelparams.resistance, v) - getgradientacceleration(params.modelparams.track, x)
    elseif mode == :Coast
        0
    elseif mode == :MaxB
        params.modelparams.umin(v)
    end
end

function solve_regular!(u0, span, p::SolverParams, seg2::Segment, shouldswitch = true, Δη = nothing)
    @unpack ρ, track, speedlimit = p.modelparams

    condition_lowspeed(u, t, int) = u[2] - 1e-2
    affect_lowspeed!(int) = terminate!(int)
    function condition_modeswitch(out, u, t, int)
        out[1] = u[3]
        out[2] = u[3] - (ρ - 1)
        out[3] = u[2] - seg2.holdspeed
        # out[3] = seg2.mode == :HoldP ? u[3] : u[3] - ρ + 1
    end
    function affect_modeswitch!(int, idx)
        if idx == 3 && int.t ≥ seg2.start && !isinf(seg2.start)
            terminate!(int)
            return
        end

        if shouldswitch
            if idx == 1 && int.p.currentmode == :Coast
                int.p.currentmode = :MaxP
                push!(switching_points, (int.t, int.p.currentmode))
            elseif idx == 2 && int.p.currentmode == :MaxB
                int.p.currentmode = :Coast
                push!(switching_points, (int.t, int.p.currentmode))
            end
        end
    end
    function affect_neg_modeswitch!(int, idx)
        if idx == 3 && int.t ≥ seg2.start && !isinf(seg2.start)
            terminate!(int)
            return
        end

        if shouldswitch
            if idx == 1 && int.p.currentmode == :MaxP
                int.p.currentmode = :Coast
                push!(switching_points, (int.t, int.p.currentmode))
            elseif idx == 2 && int.p.currentmode == :Coast
                int.p.currentmode = :MaxB
                push!(switching_points, (int.t, int.p.currentmode))
            end
        end
    end

    switching_points = []

    cb_lowspeed = ContinuousCallback(condition_lowspeed, affect_lowspeed!)

    tstops = [r[:Distance] for r in eachrow(track.waypoints)]

    cb_modeswitch = VectorContinuousCallback(condition_modeswitch, affect_modeswitch!, affect_neg_modeswitch!, 3)

    cbS = CallbackSet(cb_modeswitch, cb_lowspeed)
      
    prob = ODEProblem(odefun!, u0, span, p)
    sol = solve(prob;alg_hints = [:stiff], tstops, callback = cbS, d_discontinuities = tstops,)
    return sol, switching_points
end

function try_link(x0, seg2, initmode, modelparams::ModelParams,  across = false, vinit = modelparams.V, Δη = nothing)
    p0 = SolverParams(modelparams, initmode)
    @unpack speedlimit = modelparams

    segs = getmildsegments(modelparams)
    seg1idx = findfirst([seg.start ≤ x0 ≤ seg.finish for seg in segs])
    if isnothing(seg1idx)
        seg1 = segs[1]
    else
        seg1 = segs[seg1idx]
    end

    # Link first segment to final segment
    if across
        sol, _ = solve_regular!([0.0, vinit, 0.0], (x0, finish(modelparams.track)), p0, seg2)
        if sol.retcode == ReturnCode.Terminated
            if !isnothing(speedlimit) && sol[2,end] ≥ speedlimit(sol.t[end])
                return Inf
            end
            return -Inf
        else
            return sol[2,end] - modelparams.vf
        end
    end

    # Link first segment
    if isinf(seg2.start)
        u0 = seg2.mode == :HoldP ? [0.0, seg2.holdspeed, 0.0] : [0.0, seg2.holdspeed, modelparams.ρ - 1]
        sol, _ = solve_regular!(u0, (x0, seg2.finish), p0, seg2)    
        if sol.retcode == ReturnCode.Terminated # came to stop before finish or broke speedlimit
            if !isnothing(speedlimit) && sol[2,end] ≥ speedlimit(sol.t[end])
                return Inf
            end
            return -Inf
        elseif sol.retcode == ReturnCode.Success
            return sol[2,end] - modelparams.vᵢ
        end
    end

    # Link final segment
    if isinf(seg2.finish)
        u0 = seg1.mode == :HoldP ? [0.0, seg1.holdspeed, 0.0] : [0.0, seg1.holdspeed, modelparams.ρ - 1]
        sol, _ = solve_regular!(u0, (x0, seg2.start), p0, seg2)
        if sol.retcode == ReturnCode.Terminated # came to stop before finish or broke speedlimit
            if !isnothing(speedlimit) && sol[2,end] ≥ speedlimit(sol.t[end])
                return Inf
            end
            return -Inf
        elseif sol.retcode == ReturnCode.Success
            return sol[2,end] - modelparams.vf
        end        
    end

    # Link inner segments
    u0 = seg1.mode == :HoldP ? [0.0, seg1.holdspeed, 0.0] : [0.0, seg1.holdspeed, modelparams.ρ - 1]

    if !isnothing(Δη) # find jump in η
        sol, _ = solve_regular!(u0, (x0, seg2.finish), p0, seg2, true, Δη)
        if sol.t[end] ≈ seg2.finish
            return sign(sol[3,end]) * Inf
        elseif sol[2,end] ≤ 0.1
            return -Inf
        else
            if seg2.mode == :HoldP
                return sol[3,end]
            else # seg2.mode == :HoldR
                return sol[3,end] - (modelparams.ρ - 1)
            end
        end
    end

    sol, _ = solve_regular!(u0, (x0, seg2.finish), p0, seg2)    

    v = sol[2,:]
    η = sol[3,:]
    x = sol.t

    display(plot(sol.t, sol[2,:]))

    if !checkspeedlimit(sol, x -> speedlimit)
        return Inf
    end

    if x[end] ≈ seg2.finish
        sign(η[end]) * Inf
    elseif v[end] ≤ 0.1
        -Inf
    else
        if seg2.mode == :HoldP
            # v[end] - seg2.holdspeed
            η[end]
        elseif seg2.mode == :HoldR
            # v[end] - seg2.holdspeed
            η[end] - (modelparams.ρ - 1)
        end
    end
end

# From https://doi.org/10.1109/9.867018
trackX = [0,16e3,20e3,30e3]
trackY = [0,0,400,400]/9.81
track = HillyTrack(trackX, trackY)
myresistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)
T = 3200.0
ρ = 0.
u_max(v) = 0.125
u_min(v) = -0.25
vᵢ = 2.0
vf = 2.0

Vlim = 10.5

prob = TrainProblem(;track, resistance = myresistance, T, 
    umax = u_max, umin = u_min, ρ, vᵢ, vf, speedlimit = x -> Vlim)

chain, sol = OptimalTrainControl.solve!(prob)

plot(sol.t, sol[2,:]; color = modecolor(sol.t, chain), label = false, lw = 2)


