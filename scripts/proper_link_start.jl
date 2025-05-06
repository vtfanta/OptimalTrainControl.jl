using OptimalTrainControl
using Plots
using Roots
using StaticArrays

import OptimalTrainControl: GRAV_ACC as G

# track from Khmelnitsky, 2000
track = Track(
    40e3,
    altitude = 0.,
    x_gradient = [0., 16e3, 20e3, 24e3, 25e3, 28e3, 31e3],
    gradient = [0., (400/G)/4e3, -240/G/4e3, 0, 300/G/3e3, -180/G/3e3, 0]
)

# plot(track)

train = Train(
    v -> 0.125,
    v -> -0.25,
    (1.5e-2, 8.98e-5, 8e-5),
    0.
)

V = 12.
W = Inf
T = 2600.
eetcprob = EETCProblem(T, train, track, 2.)
simparams = EETCSimParams(eetcprob, V, W, [0.], MaxP)

port1 = Port(-Inf, 0., Coast, 12.)
port2 = Port(0., 2e3, HoldP, V)
link_sol = OptimalTrainControl.link(port1, port2, simparams)
plot(link_sol)