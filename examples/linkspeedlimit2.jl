using Plots
using OptimalTrainControl

function customlimit(x)
    if x ≤ 11e3
        10
    else
        20
    end
end

trackX = [0,10e3,14e3,25e3]
trackY = [0,0,400,400]/9.81
track = HillyTrack(trackX, trackY)
myresistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)
V = 9.625
ρ = 0.
u_max(v) = 0.125
u_min(v) = -0.25
vᵢ = 2.0
vf = 2.0

params = ModelParams(;track, resistance = myresistance, 
    V, vᵢ, vf, ρ, umax = u_max, umin = u_min, speedlimit = customlimit)

segs = getmildsegments(params)
@show segs

sol, points = linkunderspeedlimit(segs[2], segs[3], params)

plot(sol.t, sol[2,:]; color = modecolor(sol.t, points), label = false, lw = 2)


