using OptimalTrainControl
using OrdinaryDiffEq
using Plots
using Roots
using StaticArrays

track = Track(
    length = 5e3
)

train = Train(
    v -> 1/v,
    v -> -1/v,
    (1e-2, 0., 1.5e-5),
    0.6
)

T = 383.63 # seconds

prob = EETCProblem(T, train, track)
sol = solve(prob; atol=1)
plot(sol, label = false, xlabel = "Distance (m)", ylabel = "Speed (m/s)")