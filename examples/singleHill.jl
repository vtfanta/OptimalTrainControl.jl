using OptimalTrainControl
using Plots
using DifferentialEquations
using BasicInterpolators
using ForwardDiff

trackX = [0.0, 500.0, 1e3, 1.5e3]
trackY = [0.0, 0.0, 5.0, 5.0]

onehilltrack = HillyTrack(trackX, trackY)

plot!(onehilltrack, linewidth = 3)

myresistance = DavisResistance(1e-2, 0., 1.5e-5)

ρ = 0.0
V = 10.

function mycontrol(v, x)
    if v ≥ V
        resistance(myresistance, v) - getgradientacceleration(onehilltrack, x)
    else
        1. / max(v, 5.)
    end
end

"""
    findinitvalue_η(track::Track, V, initvals, maxcontrol)

Calculate the initial value of η adjoint variable by simulation from the initvals
to the speed V with applying maxcontrol(v) to achieve (v(x), η(x)) = (V, 0).
"""
function findinitvalue_η(track::Track, V, initvals, maxcontrol)
    prob = ODEProblem(function (du, u, p, t)
        v = u[2]    
    
        du[1] = inv(v)
        du[2] = (maxcontrol(v) - resistance(myresistance, v) + 
            getgradientacceleration(track, t)) * inv(v)
    end, initvals, (start(track), finish(track)))

    condition(u, t, int) = u[2] - V
    affect!(int) = terminate!(int)
    cb = ContinuousCallback(condition, affect!)

    sol = solve(prob, callback = cb)


end

prob = ODEProblem(function (du, u, p, t)
        v = u[2]
        μ₂ = u[3]

        η = μ₂ / v - 1
        ζ = η + 1 - ρ

        du[1] = inv(v)
        du[2] = (mycontrol(v, t) - resistance(myresistance, v) + 
            getgradientacceleration(onehilltrack, t)) * inv(v)
        du[3] = 1 / v^2 + u[3] * du[2] / v + 
            ForwardDiff.derivative(x -> resistance(myresistance, x), v) / v - (η > 0 ? η : 0.0) * ForwardDiff.derivative(x -> 2. / max(x, 5.), v) +
            (ζ < 0. ? -ζ : 0.) * ForwardDiff.derivative(x -> -2. / max(x, 5.), v)
    end, [0.0, 1.0, 1.0], (start(onehilltrack), finish(onehilltrack)))

sol = solve(prob; dtmax = 10)
μ₂ = sol[3,:]
v = sol[2,:]
t = sol[1,:]
x = sol.t

plot(x,v)
# plot!(x,μ₂)

# v_int = LinearInterpolator(x, v, WeakBoundaries())

# η = @. (E(myresistance, V, v_int(x)) - E(myresistance, V, V)) / (mycontrol(v_int(x), x) - resistance(myresistance, v_int(x)) + getgradientacceleration(onehilltrack, x))
# plot(η, v)