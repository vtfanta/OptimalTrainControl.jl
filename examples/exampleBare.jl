using RailDynamics
using Plots

myresistance = DavisResistance(0, 0.5, 1.3)

mymodel = AlbrechtModel(myresistance, x -> 10, x -> -10, 100)

myscenario = BasicScenario(mymodel, FlatTrack(100), 9.81)

sol = play(myscenario, [0.0, 0.5])

plot(sol[1,:], sol[2,:])