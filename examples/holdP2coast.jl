using DifferentialEquations
using NumericalIntegration
using Plots
using RailDynamics
using Roots

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
    display(plot(x,v))
    J = ψ(myresistance, V) * (Δt - Δx / V) + 
        integrate(x, resistance.(myresistance, v) .- resistance(myresistance, V))
    @show J
    J
end

function mycontrol(x, v, threshold)
    if x ≥ threshold
        2. / max(5., v)
    else
        resistance(myresistance, v) - getgradientacceleration(steephilltrack, x)
    end
end

V = 10.
localenergy(450.)


xbegin = collect(400:10:460)
Js = [localenergy(x) for x in xbegin]

scatter(xbegin, Js)