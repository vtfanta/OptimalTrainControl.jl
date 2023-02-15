using RailDynamics
using Plots

myresistance = DavisResistance(0, 0.5, 1.3)

mymodel = AlbrechtModel(myresistance, x -> 10, x -> -10, 100)

myscenario = MinimalTimeScenario(mymodel, FlatTrack(100), 9.81, [0.5], [0.5])

calculatecontrol!(myscenario)

myscenario.controllaw(3)