using BenchmarkTools
using DifferentialEquations
using NumericalIntegration
using Plots
using RailDynamics
using Roots
import ForwardDiff: derivative

trackX = [0.0, 500.0, 700., 1.5e3]
trackY = [0.0, 0.0, 5.0, 5.0]

steephilltrack = HillyTrack(trackX, trackY)

myresistance = DavisResistance(1e-2, 0., 1.5e-5)

# Calculate the local energy functional according to the formula
# from "Local energy minimization in optimal train control"
function localenergy(xpowerstart)
    prob = ODEProblem(function (du, u, p, t)
    v = u[2]

    du[1] = inv(v)
    du[2] = (mycontrol(t, v, xpowerstart) - resistance(myresistance, v) + 
        getgradientacceleration(steephilltrack, t)) * inv(v)
    end, [0.0, 10.0], (start(steephilltrack), finish(steephilltrack)))

    condition(u, t, p) = u[2] - V
    affect!(int) = terminate!(int)
    cb = ContinuousCallback(condition, affect!, nothing)

    sol = solve(prob, alg_hints = [:stiff], callback = cb)
    # display(plot(sol.t, sol[2,:]))

    # Find the end of the Power phase
    d = sol.t[end]
    a = xpowerstart
    @show a, d
    Δx = d - a
    Δt = sol(d)[1] - sol(a)[1]
    x = sol.t[sol.t .≥ a .&& sol.t .≤ d]
    v = sol[2, sol.t .≥ a .&& sol.t .≤ d]
    J = ψ(myresistance, V) * (Δt - Δx / V) + 
        integrate(x, resistance.(myresistance, v) .- resistance(myresistance, V))

    Is = [exp(-integrate(first(x,k), (ψ.(myresistance, first(v,k)) .- first(v,k).^2 .* derivative.(x->2. /max(5.,x), first(v,k))) ./ first(v,k).^3)) for k in eachindex(x)]

    decider = integrate(x, (ψ.(myresistance, v) .- ψ.(myresistance, V)) ./ v.^3 .* Is)
    @show decider
    display(plot(x,v, title = "start at $a, $decider"))

    J
end

function decider(xpowerstart)
    prob = ODEProblem(function (du, u, p, t)
    v = u[2]

    du[1] = inv(v)
    du[2] = (mycontrol(t, v, xpowerstart) - resistance(myresistance, v) + 
        getgradientacceleration(steephilltrack, t)) * inv(v)
    end, [0.0, 10.0], (start(steephilltrack), finish(steephilltrack)))

    condition(u, t, p) = u[2] - V
    affect!(int) = terminate!(int)
    cb = ContinuousCallback(condition, affect!, nothing)

    sol = solve(prob, alg_hints = [:stiff], callback = cb)
    # display(plot(sol.t, sol[2,:]))

    # Find the end of the Power phase
    d = sol.t[end]
    a = xpowerstart
    # @show a, d
    Δx = d - a
    Δt = sol(d)[1] - sol(a)[1]
    x = sol.t[sol.t .≥ a .&& sol.t .≤ d]
    v = sol[2, sol.t .≥ a .&& sol.t .≤ d]
    J = ψ(myresistance, V) * (Δt - Δx / V) + 
        integrate(x, resistance.(myresistance, v) .- resistance(myresistance, V))

    Is = [exp(-integrate(first(x,k), (ψ.(myresistance, first(v,k)) .- first(v,k).^2 .* derivative.(x->2. /max(5.,x), first(v,k))) ./ first(v,k).^3)) for k in eachindex(x)]

    decider = integrate(x, (ψ.(myresistance, v) .- ψ.(myresistance, V)) ./ v.^3 .* Is)
end

function mycontrol(x, v, threshold)
    if x ≥ threshold
        2. / max(5., v)
    else
        resistance(myresistance, v) - getgradientacceleration(steephilltrack, x)
    end
end

V = 10.
# localenergy(450.)

segs = getmildsegments(steephilltrack, V, myresistance, v -> 2. / max(5., v))

b = segs[2].finish
c = segs[3].start

## Latest possible start of MaxP
p = ModelParams(
    (u,p,t) -> 2. / max(5., u[2]),
    (u,p,t) -> resistance(myresistance, u[2]),
    (u,p,t) -> getgradientacceleration(steephilltrack, t),
    0,
    :MaxP 
)
a₋ = b
condition_hitV(u, t, int) = u[2] - V
affect_hitV!(int) = terminate!(int)
cb_hitVfrombelow = ContinuousCallback(condition_hitV, affect_hitV!, nothing)
prob_latest = ODEProblem(albrecht_odefun!, [0.0, V], (a₋, Inf), p)
d₋ = solve(prob_latest, alg_hints = [:stiff], callback = cb_hitVfrombelow).t[end] 

## Earliest possible start of MaxP
prob_earliest_proto = ODEProblem(albrecht_odefun!, [0.0, V], (c, -Inf), p)
cb_hitVfromabove = ContinuousCallback(condition_hitV, nothing, affect_hitV!)
sol_earliestproto = solve(prob_earliest_proto, alg_hints = [:stiff], callback = cb_hitVfromabove)
a₊ = sol_earliestproto.t[end]

prob_earliest_proto2 = ODEProblem(albrecht_odefun!, [0.0, 1.0], (start(steephilltrack), Inf), p)
cb_hitV = ContinuousCallback(condition_hitV, affect_hitV!)
sol_earliestproto2 = solve(prob_earliest_proto2, alg_hints = [:stiff], callback = cb_hitV)
a₊2 = sol_earliestproto2.t[end]

xbegin = collect(a₊:5:a₋)
decisions = [decider(x) for x in xbegin]
scatter(xbegin, decisions, label = "Integral condition")

switching_point = find_zero(x -> decider(x), [max(a₊, a₊2), a₋])

plot!(twinx(), steephilltrack)
vline!([switching_point], xlim = (max(a₊, a₊2), c),
    label = "Optimal switching point")
# plot!(twinx(), xbegin, Js, legend = false)