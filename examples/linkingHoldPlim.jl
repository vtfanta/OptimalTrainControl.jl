using OptimalTrainControl
using Plots

trackX = [0, 1000, 1500, 3000]
trackY = [0, 0, 10, 10]
track = HillyTrack(trackX, trackY)

limX = [500, 1200]
limY = [15, 13, 16]
# limX = [50, 120]
# limY = [13, 15, 13]
prob = TrainProblem(;track, T = 100)
setspeedlimits!(prob, limX, limY)

params = ModelParams(;track, speedlimit = prob.speedlimit, 
    speedlimitX = prob.speedlimitX,
    speedlimitY = prob.speedlimitY,
    V = 15, ρ = 0)

# params = ModelParams(; track = HillyTrack([0,16e3,20e3,24e3,25e3,28e3,31e3,40e3],[0,0,400,160,160,460,280,280]/9.81), umax = v -> 0.125, umin = v -> -0.25, resistance = DavisResistance(1.5e-2,0.127e-2/sqrt(2),0.016e-2/2), V = sqrt(2*65.43), ρ = 0.7)
# setspeedlimits!(params, [20e3, 27e3], [sqrt(2*75), sqrt(2*110), sqrt(2*70)])

segs = segmentize(params)
println(segs)

plot(params.speedlimit, start(params.track), finish(params.track); ylim = (0,30))
hline!([params.V])
plot!(twinx(), params.track; alpha = 0.5)