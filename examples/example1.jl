using RailDynamics
using Plots

myresistance = DavisResistance(0, 0.5, 1.3)

mymodel = AlbrechtModel(myresistance, 10, -10, 100)

myscenario = BasicScenario(mymodel, FlatTrack(), 9.81)

sol = play(myscenario, [0.0, 0.5])

plot(sol[1,:], sol[2,:])