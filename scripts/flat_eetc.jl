using OptimalTrainControl
using Plots
using StaticArrays

track = Track(5e3)

train = Train(
    v -> 1/v,
    v -> -1/v,
    (1e-2, 0., 1.5e-5),
    0.6
)

T = 485.63 # seconds

prob = EETCProblem(T, train, track)
sol = solve(prob; atol=2)
plot(sol, label = false, xlabel = "Distance (m)", ylabel = "Speed (m/s)")