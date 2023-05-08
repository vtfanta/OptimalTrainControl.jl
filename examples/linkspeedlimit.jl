using DifferentialEquations
import NumericalIntegration: integrate
using Plots
using OptimalTrainControl
using Roots
import ForwardDiff: derivative

function odefun!(du, u, p, x)
    t, v, η = u

    du[1] = inv(v)
    du[2] = (p.u(u, p, x) - p.r(u, p, x) + 
        p.g(u, p, x)) * inv(v)

    if p.currentmode != :MaxB
        du[3] = (ψ(myresistance, v) - v^2 * derivative(v -> p.u([t, v, η], p, x), v)) * η / v^3 + (ψ(myresistance,v) - ψ(myresistance,V))/v^3
    else
        du[3] = (ψ(myresistance, v) - v^2 * derivative(v -> p.u([t, v, η], p, x), v)) * η / v^3 + (ψ(myresistance,v) - ψ(myresistance,V))/v^3 - (1 - p.ρ) * derivative(v -> p.u([t, v, η], p, x), v) / v
    end
end

function mycontrol(u, p, x)
    if p.currentmode == :MaxP
        u_max(u[2])
    elseif p.currentmode == :Coast
        0
    elseif p.currentmode == :MaxB
        u_min(u[2])
    end
end

function solve_regular!(u0, span, p0, seg2, x0 = nothing)
    condition_lowspeed(u, t, int) = u[2] - 1e-2
    affect_lowspeed!(int) = terminate!(int)
    function condition_modeswitch(out, u, t, int)
        out[1] = u[3]
        out[2] = u[3] - (int.p.ρ - 1)
    end
    function affect_modeswitch!(int, idx)
        if idx == 1 && int.t ≥ seg2.start
            terminate!(int)
        end

        if idx == 1 && int.p.currentmode == :Coast
            int.p.currentmode = :MaxP
        elseif idx == 2 && int.p.currentmode == :MaxB
            int.p.currentmode = :Coast
        end
    end
    function affect_neg_modeswitch!(int, idx)
        if idx == 1 && int.t ≥ seg2.start
            terminate!(int)
        end

        if idx == 1 && int.p.currentmode == :MaxP
            int.p.currentmode = :Coast
        elseif idx == 2 && int.p.currentmode == :Coast
            int.p.currentmode = :MaxB
        end
    end

    function affect_speedlimit!(int)
        int.u[3] += x0
        if int.u[3] > 0
            int.p.currentmode = :MaxP
        elseif p0.ρ ≤ int.u[3] ≤ 0
            int.p.currentmode = :Coast
        else
            int.p.currentmode = :MaxB
        end
    end

    cb_modeswitch = VectorContinuousCallback(condition_modeswitch, affect_modeswitch!, affect_neg_modeswitch!, 2)
    cb_lowspeed = ContinuousCallback(condition_lowspeed, affect_lowspeed!)

    tstops = [r[:Distance] for r in eachrow(steephilltrack.waypoints)]

    if !isnothing(x0)
        cb_speedlimit = ContinuousCallback((u, t, int) -> u[2] - speedlimit, affect_speedlimit!)
        cbS = CallbackSet(cb_modeswitch, cb_lowspeed, cb_speedlimit)
    else
        cbS = CallbackSet(cb_modeswitch, cb_lowspeed)
    end
    
    prob = ODEProblem(odefun!, u0, span, p0)
    solve(prob;alg_hints = [:stiff], tstops, callback = cbS, d_discontinuities = tstops, 
        dtmax = 10)
end

function try_link(x0, seg2, initmode, linkmode = :normal, start = 0)
    if linkmode == :normal
        p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
        (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, initmode)
        sol = solve_regular!([0.0, V, 0.0], (x0, seg2.finish), p0, seg2)    

        v = sol[2,:]
        η = sol[3,:]
        x = sol.t

        @show x0
        # display(plot(x, v))

        if maximum(v) > speedlimit
            return Inf
        end

        if x[end] ≈ seg2.finish
            sign(η[end]) * Inf
        elseif v[end] ≤ 0.1
            -Inf
        else
            if seg2.mode == :HoldP
                v[end] - seg2.holdspeed
            elseif seg2.mode == :HoldR
                v[end] - seg2.holdspeed
            end
        end
    elseif linkmode == :speedlimit
        p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
        (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, initmode)
        sol = solve_regular!([0.0, V, 0.0], (start, seg2.finish), p0, seg2, x0)

        # @show sol.t[end]
        # @show sol[2,end],sol[3,end]
        # display(plot(sol.t, sol[3,:]))
        # display(vline!([seg2.start, seg2.finish]))

        if sol.t[end] ≈ seg2.finish
            sign(sol[3,end]) * Inf
        elseif sol[2,end] ≤ 0.1
            - Inf
        else
            sol[2,end] - seg2.holdspeed
        end
    end
end

function getys(γs, X)
f(xprev, xnow, yprev, γ) = tan(asin(γ / -9.81))*(xnow - xprev) + yprev
    @assert length(X) == length(γs) + 1 "Lengths don't match!"
    ret = []
    for idx in eachindex(X)
        if idx == 1
            push!(ret, 0.0)
        else
            push!(ret, f(X[idx-1],X[idx],ret[idx-1],γs[idx-1]))
        end
    end
    return ret
end

# trackX = [0,2e3,2.5e3,5e3]
# trackY = [0,0,5,5]

# steephilltrack = HillyTrack(trackX, Vector{Float64}(trackY))

# myresistance = DavisResistance(1e-2,0,1.5e-5)

# V = 25.0
# speedlimit = 25.5
# segs = getmildsegments(steephilltrack, V, myresistance, x -> 1/max(x,5))
# @show segs
# ρ = 0

trackX = [0,10e3,14e3,25e3]
trackY = [0,0,400,400]/9.81
steephilltrack = HillyTrack(trackX, trackY)
myresistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)
V = 9.625
ρ = 0.
u_max(v) = 0.125
u_min(v) = -0.25
vᵢ = 2.0
vf = 2.0

speedlimit = 10.2

segs = getmildsegments(steephilltrack, V, myresistance, u_max)
@show segs

startingmode = :MaxP
targetseg = segs[3]

xopt = find_zero(x -> try_link(x, targetseg, startingmode), (segs[2].start, segs[2].finish-1))
p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, startingmode)
sol = solve_regular!([0.0,V,0.0], (xopt, targetseg.finish), p0, targetseg)
if maximum(sol[2,:]) ≈ speedlimit
    println("CONDITION!")
    Δη = find_zero(η -> try_link(η, targetseg, startingmode, :speedlimit, xopt), (0, 5))
    p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, startingmode)
    sol = solve_regular!([0.0,V,0.0], (xopt, targetseg.finish), p0, targetseg, Δη)
end

plot(sol.t, sol[2,:]; color = [e ≥ 0 ? :green : :grey for e in sol[3,:]], lw = 3, label = false)
plot!(twinx(), steephilltrack; alpha = 0.5, label = false)
hline!([speedlimit]; label = false)