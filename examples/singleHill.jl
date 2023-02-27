using RailDynamics
using Plots
using DifferentialEquations

trackX = [0.0, 500.0, 1e3, 1.5e3]
trackY = [0.0, 0.0, 30.0, 30.0]

onehilltrack = HillyTrack(trackX, trackY)

getgrade(onehilltrack, 802.0)

plot!(onehilltrack, linewidth = 3)

myresistance = DavisResistance(1e-2, 0., 1.5e-5)

function mycontrol(v)
    return 1 / max(v, 5)
end

prob = ODEProblem(function (du, u, p, t)
        du[1] = inv(u[2])
        du[2] = (mycontrol(u[2]) - resistance(myresistance, t) + 
            getgradientacceleration(onehilltrack, t)) * inv(u[2])
    end, [0.0, 1.0], (start(onehilltrack), finish(onehilltrack)))



