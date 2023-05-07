using LaTeXStrings
using OptimalTrainControl
using Plots

trackX = [0,2e3,3e3,5e3]
trackY = [0,0,60,65]

track = HillyTrack(trackX, trackY)

params = ModelParams(track = track, V = 9, vᵢ = 12.1)

chain, sol = OptimalTrainControl.link(segs[1], segs[2],params)
# prob = TrainProblem(;track, T = 700, vᵢ = 9.71, vf = 12)
# chain, sol = solve!(prob)
# plot(sol.t, sol[2,:]; color = modecolor(sol.t, chain))