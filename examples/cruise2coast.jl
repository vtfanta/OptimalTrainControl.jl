using DifferentialEquations
using OptimalTrainControl
using Roots

function decider(xcruisestart)
    prob = ODEProblem(function (du, u, p, t)
    v = u[2]

    du[1] = inv(v)
    du[2] = ( - resistance(myresistance, v) + 
        getgradientacceleration(mytrack, t)) * inv(v)
    end, [0.0, 10.0], (xcruisestart, finish(mytrack)))

    condition(u, t, p) = u[2] - V
    affect!(int) = terminate!(int)
    cb = ContinuousCallback(condition, nothing, affect!)

    sol = solve(prob, alg_hints = [:stiff], callback = cb)
    # display(plot(sol.t, sol[2,:]))

    # Find the end of the Power phase
    d = sol.t[end]
    a = xcruisestart
    # @show a, d
    Δx = d - a
    Δt = sol(d)[1] - sol(a)[1]
    x = sol.t[sol.t .≥ a .&& sol.t .≤ d]
    v = sol[2, sol.t .≥ a .&& sol.t .≤ d]
    J = ψ(myresistance, V) * (Δt - Δx / V) + 
        integrate(x, resistance.(myresistance, v) .- resistance(myresistance, V))

    Is = [exp(-integrate(first(x,k), (ψ.(myresistance, first(v,k))) ./ first(v,k).^3)) for k in eachindex(x)]

    decider = integrate(x, (ψ.(myresistance, v) .- ψ.(myresistance, V)) ./ v.^3 .* Is)
end

trackX = [0.0, 500.0, 700., 1.5e3]
trackY = [0.0, 0.0, -1.0, -1.0]
mytrack = HillyTrack(trackX, trackY)

myresistance = DavisResistance(1e-2, 0., 1.5e-4)

V = 10.0

b = 500
c = 700

## Latest possible start to coast phase
p = ModelParams(
    (u,p,t) -> 0,
    (u,p,t) -> resistance(myresistance, u[2]),
    (u,p,t) -> getgradientacceleration(mytrack, t),
    0,
    :MaxP
)
prob_latest = ODEProblem(albrecht_odefun!, [0.0, V], (b, Inf), p)
condition_hitV(u, t, int) = u[2] - V
affect_hitV!(int) = terminate!(int)
cb_hitV = ContinuousCallback(condition_hitV, affect_hitV!)
sol_latest = solve(prob_latest, alg_hints = [:stiff], callback = cb_hitV)
plot(sol_latest.t, sol_latest[2,:])
d₊ = sol_latest.t[end]
a₊ = b

## Earliest possible start
prob_earliest = ODEProblem(albrecht_odefun!, [0.0, V], (c, -Inf), p)
sol_earliest = solve(prob_earliest, alg_hints = [:stiff], callback = cb_hitV)
plot(sol_earliest.t, sol_earliest[2,:])
##TODO it is possible this will fail and the velocity would go below 0!
a₋ = sol_earliest.t[end]

xbegin = collect(a₋:5:a₊)
deciders = decider.(xbegin)
scatter(xbegin, deciders, label = "Integral condition")

switching_point = find_zero(x -> decider(x), [a₋, a₊])

plot!(twinx(), mytrack)

vline!([switching_point], label = "Optimal switching point",
    xlim = (a₋, c))