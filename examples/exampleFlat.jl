using RailDynamics
using Plots
using NumericalIntegration

myresistance = DavisResistance(1e-2, 0., 1.5e-5)

mymodel = AlbrechtModel(myresistance, v -> inv(max(5., v)), v -> -inv(max(5., v)), 1e3)

myscenario = MinimalTimeScenario(mymodel, FlatTrack(30e3), 9.81, [0.5], [0.5])

ts, bs = calculatecontrol!(myscenario);

sol = play(myscenario)

plot(sol, label = "Result")
plot!(ts, label = "Throttle")
plot!(reverse(bs.t), reverse(bs[1,:]), label = "Brake")

# Total time calculation
velocity = sol[1,:]
T = integrate(sol.t, inv.(velocity))
t = cumul_integrate(sol.t, inv.(velocity))