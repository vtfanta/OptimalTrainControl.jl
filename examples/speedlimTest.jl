using OptimalTrainControl
using Plots

trackX = [0,2e3,2.5e3,5e3]
trackY = [0,0,5,5]

track = HillyTrack(trackX, trackY)

res = DavisResistance(1e-2,0,1.5e-5)
V = 25.0
speedlimit = 25.5

prob = TrainProblem(;track, T = 1300, resistance = res, umax = v -> 1 / max(5,v))
setspeedlimits!(prob, [],[speedlimit])

params = ModelParams(prob.umax, prob.umin, prob.resistance, prob.ρ, prob.track, V, prob.vᵢ, prob.vf, prob.speedlimit)

segs = getmildsegments(params)

xopt = 1355.5753692062572

l = OptimalTrainControl.link(segs[2], segs[3], params)
plot(l[1].t,l[1][2,:])

# solve!(prob)
# chain, sol = solve!(prob)

# plot(sol.t, sol[2,:]; color = modecolor(sol.t, chain),
    # lw = 3, label = false)