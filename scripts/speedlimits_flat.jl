using OptimalTrainControl
using Plots

train = Train(
    v -> 1/v,
    v -> -1/v,
    (1e-2, 0., 1.5e-5),
    0.6
)

track = Track(
    length = 5e3,
    altitude = 10.,
    x_speedlimit = [0., 1e3, 2.5e3],
    speedlimit = [9., 30., 8.]
)

segments = hold_segments!(EETCProblem(500., train, track), 10., 12.)
plot(segments)