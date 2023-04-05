using DifferentialEquations
import NumericalIntegration: integrate
using Plots
using RailDynamics
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
    p0 = ModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, initmode)
    sol = solve_regular!([0.0, V, 0.0], (x0, seg2.finish), p0, seg2)    

    v = sol[2,:]
    η = sol[3,:]
    x = sol.t

    @show x[end]

    if x[end] ≈ seg2.finish
        sign(η[end]) * Inf
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

trackX = [-1e3,-0.8e3,-0.2e3,0,1e3,1.8e3,1.9e3,3.1e3,4.2e3,4.6e3,5e3]
trackY = getys([0,0.03,0,0,-0.03,0,0.03,0,-0.04,0], trackX)

steephilltrack = HillyTrack(trackX, Vector{Float64}(trackY))

myresistance = DavisResistance(1e-2,0,1.5e-5)

V = 25.0
segs = getmildsegments(steephilltrack, V, myresistance, x -> 1/max(x,5))
@show segs
ρ = 0

startingmode = :Coast
targetseg = segs[3]

xopt = find_zero(x -> try_link(x, targetseg, startingmode), (segs[2].start, segs[2].finish-1))
p0 = ModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, startingmode)
sol1 = solve_regular!([0.0,V,0.0], (xopt, targetseg.finish), p0, targetseg)
plot(sol1.t, sol1[2,:]; color = [e ≥ 0 ? :green : :grey for e in sol1[3,:]], lw = 3, label = false,
    ylabel = "Speed (m/s)")

startingmode = :MaxP
targetseg = segs[5]

xopt = find_zero(x -> try_link(x, targetseg, startingmode), (segs[3].start, segs[3].finish-1))
p0 = ModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, startingmode)
sol2 = solve_regular!([0.0,V,0.0], (xopt, targetseg.finish), p0, targetseg)
plot!(sol2.t, sol2[2,:]; color = [e ≥ 0 ? :green : :grey for e in sol2[3,:]], lw = 3, label = false)

startingmode = :MaxP
targetseg = segs[6]

xopt = find_zero(x -> try_link(x, targetseg, startingmode), (segs[5].start, segs[5].finish-1))
p0 = ModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, startingmode)
sol3 = solve_regular!([0.0,V,0.0], (xopt, targetseg.finish), p0, targetseg)
plot!(sol3.t, sol3[2,:]; color = [e ≥ 0 ? :green : :grey for e in sol3[3,:]], lw = 3, label = false)

plot!(twinx(), steephilltrack; alpha = 0.5, size = (1600/2,900/2))

plot!([sol1.t[end],sol2.t[1]],[V,V]; color = :blue, lw = 3, label = false)
plot!([sol2.t[end],sol3.t[1]],[V,V]; color = :blue, lw = 3, label = false)