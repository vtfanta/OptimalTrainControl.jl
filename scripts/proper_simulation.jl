using DiffEqCallbacks
using OptimalTrainControl
using OrdinaryDiffEq
using Plots
using StaticArrays

# use this script to finetune the process of running simulation of regular phases (MaxP, Coast, MaxB)
# - ensure proper transition between phases
# - ensure proper calculation and saving of η
# - ensure proper calculations of Es (necessary for η)

# Test on examples in 10.1016/j.automatica.2013.07.008 and 10.1016/j.trb.2015.07.024

# define the track
γ2grade(γ) = (-γ / 9.81) / sqrt(1 - (γ / 9.81)^2)
γs = [0, -0.2, 0, -0.26, 0]
# track = Track(
#     length = 7.2e3,
#     altitude = 0.,
#     x_gradient = [0., 5e3, 5.6e3, 5.8e3, 6.3e3], 
#     gradient = γ2grade.(γs),
# )

# track = Track(
#     length = 9.5e3,
#     x_gradient = [0., 5e3, 5.6e3, 5.8e3, 6.4e3],
#     gradient = γ2grade.([0., -0.2, 0., 0.2, 0.])
# )

# track = Track(
#     length = 7e3,
#     x_gradient = [0., 1.5e3, 4e3, 4.4e3, 5.5e3],
#     gradient = γ2grade.([0., -0.025, 0., -0.03, 0.])
# )
track = OptimalTrainControl.Track(
    4e3,
    x_gradient = [0., 1e3, 1.8e3, 1.9e3, 3.1e3],
    gradient = γ2grade.([0., -0.03, 0., 0.03, 0.])
)

# define the train
# train = Train(
#     v -> 2/v,
#     v -> -3/v,
#     (6.75e-3, 0., 5e-5),
#     0.0
# )

# V = 20.
V = 25.

train = Train(
    v -> 1. / max(v, 5.),
    v -> - 3. / max(v, 5),
    (1e-2, 0., 1.5e-5),
    0.
)

# setup simulation
s0 = SA[0., V]
# xspan = (4662, 6664)
xspan = (837., 3.5e3)

# for debugging only
function _get_initial_E(η₀, params::EETCSimParams) 
    phase = params.current_phase
    V = params.V
    W = params.W
    train = params.eetcprob.train
    track = params.eetcprob.track
    v0 = params.eetcprob.initial_speed
    if phase == MaxP
        Es1 = η₀ * (train.U̅(v0) - r(train, v0) + g(track, track.x_gradient[1])) - E(train, V, v0)
    elseif phase == Coast
        Es1 = η₀ * (-r(train, v0) + g(track, track.x_gradient[1])) - E(train, V, v0)
    elseif phase == MaxB
        Es1 = η₀ * (train.U̲(v0) - r(train, v0) + g(track, track.x_gradient[1])) - train.ρ * E(train, W, v0)
    end
    Es1
end

# define the EETC problem
total_time = 500.
prob = EETCProblem(total_time, train, track)
prob.initial_speed = V

simparams = EETCSimParams(prob, V, OptimalTrainControl.calculate_W(prob, V), [0.], MaxP)
# guess initial Es[1]; during solution of EETC problem, this needs to be optimized for
Es1 = _get_initial_E(0., simparams)
simparams = EETCSimParams(prob, V, OptimalTrainControl.calculate_W(prob, V), [Es1], MaxP)

# simulate
otc_sol = simulate_regular_forward(simparams, xspan, s0)
odesol = otc_sol.odesol
η = otc_sol.η

# plot(odesol.t, modes)
plot(odesol.t, η)
# plot(η, odesol[2,:])

## Profiling
# function prof_g(n, simparams, xspan, s0)
#     for _ in 1:n
#         simulate_regular_forward(simparams, xspan, s0)
#     end
# end

# @profview_allocs prof_g(100, simparams, xspan, s0)