using OptimalTrainControl
using Plots

track = Track(
    5e3,
    altitude = 10.,
    x_gradient = [0., 2e3, 3e3, 4e3],
    gradient = [0., 35/1e3, 0., -35/1e3],
    x_speedlimit = [0., 1e3, 2.5e3],
    speedlimit = [9., 30., 15.]
)

train = Train(
    v -> 1/v,
    v -> -1/v,
    (1e-2, 0., 1.5e-5),
    0.6
)

prob = EETCProblem(500., train, track)

ports = hold_segments!(prob, 10.)

plot(ports)
ylims!(0., 12)
plot!(twinx(), track)