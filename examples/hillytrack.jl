using OptimalTrainControl
using Plots

# From https://doi.org/10.1109/9.867018
trackX = [0,16e3,20e3,24e3,25e3,28e3,31e3,40e3]
trackY = [0,0,400,160,160,460,280,280]/9.81
track = HillyTrack(trackX, trackY)
myresistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)

T = 2600.0
ρ = 0.
u_max(v) = 0.125
u_min(v) = -0.25
vᵢ = 2.0
vf = 2.0

prob = TrainProblem(;track = subtrack(track, 11e3), resistance = myresistance, T, 
    umax = u_max, umin = u_min, ρ, vᵢ, vf)
points, sol = solve!(prob)
plot(sol.t, sol[2,:]; color = modecolor(sol.t, chain), label = false, lw = 3, ylabel = "Speed (m/s)")
plot!(twinx(), track; alpha = 0.5, xlabel = "")