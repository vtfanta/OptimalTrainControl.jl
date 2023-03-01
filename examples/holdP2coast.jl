using RailDynamics
using Roots
using Plots
using DifferentialEquations
using NumericalIntegration

trackX = [0.0, 500.0, 700., 1.5e3]
trackY = [0.0, 0.0, 5.0, 5.0]

steephilltrack = HillyTrack(trackX, trackY)

myresistance = DavisResistance(1e-2, 0., 1.5e-5)

# Calculate the local energy functional according to the formula
# from 
function localenergy(xpowerstart)

end

function mycontrol(x, v, threshold)
    if x ≥ threshold
        2. / max(5., v)
    else
        resistance(myresistance, v) - getgradientacceleration(steephilltrack, x)
    end
end

threshold = 450.

prob = ODEProblem(function (du, u, p, t)
        v = u[2]

        du[1] = inv(v)
        du[2] = (mycontrol(t, v, threshold) - resistance(myresistance, v) + 
            getgradientacceleration(onehilltrack, t)) * inv(v)
    end, [0.0, 10.0], (start(steephilltrack), finish(steephilltrack)))

