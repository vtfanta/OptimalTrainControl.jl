using DifferentialEquations
using RailDynamics

prob = ODEProblem(function (du, u, p, t)
        du[1] = 1.1u[1]
        du[2] = 0 end, [1.0,-1.0], (0.0, 1.0))

