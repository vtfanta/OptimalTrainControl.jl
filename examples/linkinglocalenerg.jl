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
        1. / max(u[2],5)
    elseif p.currentmode == :Coast
        0
    elseif p.currentmode == :MaxB
        - 3. / u[2]
    end
end

function solve_regular!(u0, span, p0, seg2)
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
    
    prob = ODEProblem(odefun!, u0, span, p0)
    solve(prob;alg_hints = [:stiff], tstops, callback = cbS, d_discontinuities = tstops, 
        dtmax = 10)
end

function try_link(x0, seg2, initmode)
    p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, initmode)
    sol = solve_regular!([0.0, V, 0.0], (x0, seg2.finish), p0, seg2)    

    v = sol[2,:]
    η = sol[3,:]
    x = sol.t

    # @show x[end]

    if x[end] ≈ seg2.finish
        sign(η[end])*Inf
    elseif v[end] ≤ 0.1
        -Inf
    else
        v[end] - V
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

function local_energy(u0, x0, p0, seg2)
    sol = solve_regular!(u0, (x0, seg2.finish), p0, seg2)
    v = sol[2,:]
    x = sol.t
    J = integrate(x, E.(myresistance, V, v) .- E(myresistance, V, V))
end

trackX = [0,1.5e3,4e3,4.4e3,5.5e3,7e3]
trackY = getys([0.0,-0.025,0,-0.03,0], trackX)

steephilltrack = HillyTrack(trackX, Vector{Float64}(trackY))

myresistance = DavisResistance(1e-2, 0., 1.5e-5)

V = 25.0
segs = getmildsegments(steephilltrack, V, myresistance, x -> 1/max(x,5))
@show segs
ρ = 0

initmode = :MaxP
targetseg = segs[4]
xopt = find_zero(x -> try_link(x, targetseg, initmode), (segs[2].start, segs[2].finish-1))

p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, initmode)
sol = solve_regular!([0.0, V, 0.0], (xopt, targetseg.finish), p0, targetseg)
plot(sol.t, sol[2,:])
hline!([V])
plot!(twinx(), steephilltrack; alpha = 0.5)

# p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
#     (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, initmode)
# xs = collect(segs[2].start:10:segs[2].finish-1)
# vals = [local_energy([0.0,V,0.0], x, p0, targetseg) for x in xs]
# scatter(xs, vals)