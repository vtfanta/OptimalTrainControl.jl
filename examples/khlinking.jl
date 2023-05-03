using DifferentialEquations
using Plots
using OptimalTrainControl
using Roots
import ForwardDiff: derivative

function odefun!(du, u, p, x)
    K, Ψ = u

    atr = (Ψ ≥ 1 ? Ψ - 1 : 0)
    ab = (Ψ ≥ α ? 0 : α - Ψ)

    du[1] = mycontrol(u, p, x) - w(K) - derivative(P, x)
    du[2] = (Ψ * R(K) - R(Kₛ)) / K^(3/2) - atr * derivative(gtr, K) - ab * derivative(gb, K)
end

function mycontrol(u, p, x)
    K = u[1]
    if p.currentmode == :MaxP
        gtr(K)
    elseif p.currentmode == :Coast
        0
    elseif p.currentmode == :MaxB
        -gb(K)
    end
end

function solve_regular!(u0, xbegin, p0, seg2)
    condition_lowspeed(u, t, int) = u[1] - 5.0
    affect_lowspeed!(int) = terminate!(int)
    function condition_modeswitch(out, u, t, int)
        out[1] = u[2] - 1.0
        out[2] = u[2] - α
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
        # @show int.t, int.p.currentmode
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
        # @show int.t, int.p.currentmode
    end

    cb_modeswitch = VectorContinuousCallback(condition_modeswitch, affect_modeswitch!, affect_neg_modeswitch!, 2)
    cb_lowspeed = ContinuousCallback(condition_lowspeed, affect_lowspeed!)

    tstops = [r[:Distance] for r in eachrow(steephilltrack.waypoints)]
    
    cbS = CallbackSet(cb_modeswitch, cb_lowspeed)
    
    prob = ODEProblem(odefun!, u0, (xbegin, seg2.finish), p0)
    solve(prob;alg_hints = [:stiff], tstops, callback = cbS, d_discontinuities = tstops, 
        dtmax = 10)
end

function try_link(x0, seg2, initmode)
    p0 = ModelParams(nothing, nothing, nothing, nothing, initmode)
    sol = solve_regular!([Kₛ, 1.0], x0, p0, seg2)    

    K = sol[1,:]
    Ψ = sol[2,:]
    x = sol.t

    # @show x[end]

    if x[end] ≈ seg2.finish
        Inf
    elseif x[end] ≤ seg2.start
        -Inf
    else
        K[end] - Kₛ
    end
end

f(xprev, xnow, yprev, γ) = tan(asin(γ / -9.81))*(xnow - xprev) + yprev

function getys(γs, X)
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

w(K) = 1.5e-2 + 0.127e-2√K + 0.016e-2K
R(K) = K^(3/2)*derivative(w, K)

gtr(K) = 0.125
gb(K) = 0.25

Kₛ = 63.27
α = 0
myresistance = DavisResistance(1.5e-2, 0.127e-2/√2, 0.016e-2/2)

trackX = [0,16e3,20e3,24e3,25e3,28e3,31e3,40e3]
trackY = [0,0,400,160,160,460,280,280]
steephilltrack = HillyTrack(trackX, trackY)

segs = getmildsegments(steephilltrack, sqrt(2Kₛ), myresistance, x -> gtr(x))

P(s) = steephilltrack.interpolator(s)

x0 = segs[2].finish-1000
targetseg = segs[4]
startphase = :MaxP

p0 = ModelParams(nothing, nothing, nothing, nothing, startphase)
sol = solve_regular!([Kₛ, 1.0], 15757.489720241621, p0, targetseg)

try_link(x0, targetseg, startphase)

plot(sol.t, sol[1,:])
plot!(twinx(), steephilltrack; alpha = 0.5)

# xs = x0:10:segs[2].finish
# vals = clamp.(try_link.(xs, targetseg, startphase), -Kₛ, Kₛ)
# scatter(xs, vals)
