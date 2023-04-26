using RailDynamics
using Plots

# From https://doi.org/10.1109/9.867018
trackX = [0,16e3,20e3,24e3,25e3,28e3,31e3,40e3]
trackY = [0,0,400,160,160,460,280,280]/9.81
track = HillyTrack(trackX, trackY)
myresistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)
T = 3621.0
ρ = 0
u_max(v) = 0.125
u_min(v) = -0.25
vᵢ = 2.0
vf = 2.0

prob = TrainProblem(;track, resistance = myresistance, T = T, 
    umax = u_max, umin = u_min, ρ, vᵢ, vf)
solve!(prob)
chain, sol = prob.switchingpoints, prob.states 

display(plot(sol.t, sol[2,:]; color = modecolor(sol.t, chain), label = false, lw = 3))
@show sol[1,end] - T

