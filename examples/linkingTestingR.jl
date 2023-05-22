using OptimalTrainControl
using Plots

# trackX = [0,16e3,20e3,24e3,25e3,28e3,31e3,40e3]
# trackY = [0,0,400,160,160,460,280,280]/9.81
trackX = [0,15e3,24e3,35e3]
trackY = [0,0,-65,-65]

track = HillyTrack(trackX, trackY)
myresistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)
V = 15.5
ρ = 0.8
u_max(v) = 0.125
u_min(v) = -0.25
vᵢ = 15.0
vf = 16.0

# modelparams = OptimalTrainControl.ModelParams(u_max, u_min, myresistance, ρ, track, V, vᵢ, vf)
# segs = OptimalTrainControl.getmildsegments(modelparams)
# l = OptimalTrainControl.link(segs[3],segs[4],modelparams)
# println(l[2])
# display(plot(l[1].t,l[1][2,:]; color = modecolor(l[1].t, l[2]), lw = 2, label = false))
# display(plot!(twinx(), track; alpha = 0.5))

# chain, sol = OptimalTrainControl.findchain(segs, modelparams)
# plot(sol.t, sol[2,:]; color = modecolor(sol.t, chain))

prob = TrainProblem(track = track, T = 2600, ρ = ρ,
    vᵢ = vᵢ, vf = vf)
solve!(prob)

chain, sol = prob.switchingpoints, prob.states 
plot(sol.t, sol[2,:]; color = modecolor(sol.t, chain), lw = 2, label = false,
ylabel = "Speed (m/s)", size=(750,300), margin=5mm, ylim=(1,22), xlabel="Distance (m)")
plot!(twinx(), track; alpha = 0.5,xlabel = "")

savefig("examples/figures/resultHoldRandSpeeds.pdf")