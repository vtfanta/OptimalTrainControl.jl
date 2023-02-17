using RailDynamics
using Plots

myresistance = DavisResistance(1e-2, 0., 1.5e-5)

mymodel = AlbrechtModel(myresistance, v -> inv(max(5., v)), v -> -inv(max(5., v)), 1e3)

myscenario = OptimalScenario(mymodel, FlatTrack(30e3), 9.81, [0.0, 0.5], [2000., 0.5])

