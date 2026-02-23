using OptimalTrainControl
using Plots

# flat track with constant speed limit
track = Track(
    10e3;  # 10 km long
    altitude = 0.0,
    x_gradient = [0.0],
    gradient = [0.0],
    x_speedlimit = [0.0, 4e3],
    speedlimit = [18.0, 20.0]  # 18 m/s speed limit
)

P = 2.
a, b, c = 6.75e-3, 0., 5e-5
train = Train(
    v -> P / v,
    v -> -P / v,
    (a, b, c),
    0.
)

V = 12.0 # cruising speed
prob = EETCProblem(8*60., train, track, 1.)

ports = hold_segments!(prob, V)
start_port = Port(-Inf, 0.0, MaxP, 1.)

link_cond = OptimalTrainControl.LinkPortConditionModule.makeLinkPortCondition(ports[1], prob)
W = OptimalTrainControl.calculate_W(prob, V)

OptimalTrainControl.root_f_start(-2., prob, V, W, link_cond)
# OptimalTrainControl.link_start(ports[1], EETCSimParams(prob, V, MaxP))