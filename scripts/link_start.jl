using OptimalTrainControl
using OrdinaryDiffEq
using DiffEqCallbacks
using Plots
using Roots
using StaticArrays

include("./eta_utils.jl")

train = Train(
    v -> 3/v,
    v -> -3/v,
    (6.75e-3, 0., 5e-5),
    0.
)

track = Track(;
    length = 7.2e3,
    x_gradient = [0., 5e3, 5.6e3, 5.8e3, 6.3e3],
    gradient = asin.([0, -0.2, 0, -0.26, 0] ./ -9.81)
)

# want to start from 1 m/s and reach cruising speed
xspan = SA[0., 1e3]
v0 = 1.
cruising_speed = 12.
target_speed = 11.

# function f(E0)
E0 = -0.0285
eetcprob = EETCProblem(;
    train,
    track = track,
    current_phase = MaxP,   # tell solver to start with MaxP; later this should be inferred based on initial speed
    initial_speed = v0,
    T = 1e3,
    Es = [E0]
)

# I want to search for such initial η (Es[1]), that the train reaches cruising_speed when η=0.
# For this example, correct is -0.2835 (if target_speed == cruising_speed)
# starting port
start_port = Port(;
    start = -Inf,
    finish = 0.,
    mode = HoldP,
    speed = v0)

odeprob = ODEProblem(_odefun, SA[0., v0], xspan, eetcprob;
    tstops = track.x_gradient)

# saving η
saved_vals = SavedValues(Float64, Float64)  # time type, savedval type
saving_cb = SavingCallback((s,x,int) -> calculate_η(s,x,int.p,cruising_speed), saved_vals, save_everystep = true)

# stop when reaching cruising speed
targetspeed_cond(s, x, int) = s[2] - target_speed
targetspeed_aff(int) = terminate!(int)
targetspeed_cb = ContinuousCallback(targetspeed_cond, targetspeed_aff, save_positions=(false, false))
# save_positions saves before and after callback application if set to (true, true)
odesol = OrdinaryDiffEq.solve(odeprob, AutoVern7(Rodas5());
    callback = CallbackSet(saving_cb, targetspeed_cb),
    d_discontinuities = track.x_gradient,
    tstops = track.x_gradient,
    dtmax = 3.)

# vector of η values, x points are the same as odesol.t
η = saved_vals.saveval
η[end]
# end

# z = find_zero(f, 0)