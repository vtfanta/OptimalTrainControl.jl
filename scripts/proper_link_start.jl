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

function root_f(E1, eetcprob, V, W)
# E1 = 0.
    Es = [E1]
    simparams = EETCSimParams(eetcprob, V, W, Es, MaxP)
    otc_sol = simulate_regular_forward(simparams, (0., 1e3), SA[0., eetcprob.initial_speed]);
    idx = searchsortedfirst(otc_sol.odesol[2,:], V)
    if idx > length(otc_sol.odesol[2,:])
        -Inf
    else
        otc_sol.η[idx]
    end
end

e = Roots.find_zero(e -> root_f(e, eetcprob, V, W), [-30.,30.], atol=1e-6)
# EETCSimParams(eetcprob, V, W, [-0.1], MaxP)
# otc_sol = simulate_regular_forward(simparams, (0., 1e3), SA[0., eetcprob.initial_speed]);
# plot(otc_sol)

##
simparams = EETCSimParams(eetcprob, V, W, [e], MaxP)
con = OptimalTrainControl.LinkPortConditionModule.LinkPortCondition(0., 3e3, V)
link_sol = OptimalTrainControl.simulate_link_forward(0., con, simparams, SA[0., eetcprob.initial_speed])

using Plots
plot(link_sol)