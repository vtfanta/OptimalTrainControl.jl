using OptimalTrainControl
using Plots

train = Train(
    U̅ =  v -> 3/v,
    U̲ =  v -> -3/v,
    r = (6.75e-3, 0., 5e-5)
)

track = Track(
    length = 3e3,
    altitude = 100.,
    x_gradient = [0.0, 1e3, 1.7e3],
    gradient = [2e-3, 0., 1e-3]
)

prob = TOTCProblem(train, track)

sol = solve(prob)

plot(sol, label = false, ylabel = "Speed (m/s)", xlabel = "Distance (m)")
plot!(twinx(), track, ylabel = "Altitude (m)")

# plot(sol.control, 0, track.length, label = false, ylabel = "Control", xlabel = "Distance (m)")