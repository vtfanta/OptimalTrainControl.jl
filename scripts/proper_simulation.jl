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

const VECTOR_CB_LENGTH = 3

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

"""
In-place version of the right-hand side of the ODE system for the EETC problem.

Params are in a EETCSimParams struct.
"""
function _odefun(s::A, p::EETCSimParams, x) where {T<:Real, A<:AbstractArray{T,1}}
    _, v = s

    eetcprob, phase = p.eetcprob, p.current_phase
    if phase == MaxP
        u = eetcprob.train.U̅(v)
    elseif phase == Coast
        u = 0.
    elseif phase == MaxB
        u = eetcprob.train.U̲(v)
    end
    # @show u

    ds1 = 1/v
    ds2 = (u - OptimalTrainControl.r(eetcprob.train, v) + OptimalTrainControl.g(eetcprob.track, x)) / v
    SA[ds1, ds2]
end 

function calculate_η(s, x, p::EETCSimParams{T}) where {T<:Real}
    v = s[2]

    eetcprob, V, W, Es, phase = p.eetcprob, p.V, p.W, p.Es, p.current_phase
    train = eetcprob.train
    track = eetcprob.track
    
    # if trying to evaluate η at invalid point (Es not updated yet), return same sign to prevent root
    x_index = searchsortedfirst(track.x_gradient, x)
    if x_index > length(track.x_gradient) || track.x_gradient[x_index] != x
        if length(p.Es) < x_index - 1   # Es invalid, return sign of last valid value
            x = track.x_gradient[x_index - 1] - 1e-2
        else
            # valid Es
        end
    else
        if length(p.Es) < x_index
            x = track.x_gradient[x_index] - 1e-2
        else
            # valid Es
        end
    end

    val::T = zero(T)
    if phase == MaxP
            val = (E(train, V, v) + Es[end]) / 
            (train.U̅(v) - OptimalTrainControl.r(train, v) + g(track, x))
    elseif phase == Coast
        val = (E(train, V, v) + Es[end]) / 
            (- OptimalTrainControl.r(train, v) + g(track, x))
    elseif phase == MaxB
        # have to calculate ζ, the Fs are also in the EETC.Es array
        # TODO can E be found such that I dont have to store F instead?
        val = (train.ρ*E(train, W, v) + Es[end]) /
            (train.U̲(v) - OptimalTrainControl.r(train, v) + g(track, x))
        # η = ζ + ρ - 1
        val += train.ρ - 1.
    end
    return val
end

function simulate_regular_forward(sparams::EETCSimParams, xspan::Tuple{T, T}, s0::SVector{2, T}) where {T<:Real}
    
    # to prevent mutating
    simparams = deepcopy(sparams)

    function mode_switch_condition(out, s, x, integrator)
        _, v = s
        p = integrator.p
        η = calculate_η(s, x, p)

        # @show x, integrator.t, integrator.tprev
        # @show any(x .≈ p.eetcprob.track.x_gradient)
        # @show calculate_η(s, x, p), calculate_η(s, integrator.t, p), calculate_η(s, integrator.tprev, p)
        out[1] = η  # boundary between MaxP and Coast
        out[2] = η - (p.eetcprob.train.ρ - 1.)    # boundary between Coast and MaxB
        # TODO should I add terminating condition?
        out[3] = v - 0.1    # terminate when speed is too low
    end

    function mode_switch_affect_pos!(integrator, index) # gets applied when crossing 0 from negative to positive
        if index == 1   # Coast to MaxP
            integrator.p.current_phase = MaxP
        elseif index == 2   # MaxB to Coast
            integrator.p.current_phase = Coast
            # append to Es since switching from ζ to η
            v = integrator.u[2]
            Es = integrator.p
            eetcprob = integrator.p.eetcprob
            train = eetcprob.train
            track = eetcprob.track
            η₀ = train.ρ - 1.  # know this from analytic derivation
            
            newE = η₀ * (-r(train, v) + g(track, integrator.t)) - E(train, integrator.p.V, v)
            push!(Es, newE)
        end
    end

    function mode_switch_affect_neg!(integrator, index) # gets applied when crossing 0 from positive to negative
        if index == 1   # MaxP to Coast
            integrator.p.current_phase = Coast
        elseif index == 2   # Coast to MaxB 
            v = integrator.u[2]
            Es = integrator.p.Es
            eetcprob = integrator.p.eetcprob
            train = eetcprob.train
            track = eetcprob.track
            η₀ = train.ρ - 1.  # know this from analytic derivation

            integrator.p.current_phase = MaxB

            # append to Es since switching from η to ζ (it's actually F)
            newF = (η₀ - train.ρ + 1.) * (train.U̲(v) - r(train, v) + g(track, integrator.t)) - train.ρ * E(train, integrator.p.W, v)
            push!(Es, newF)
        elseif index == 3   # terminate because of low speed
            terminate!(integrator)
        end
    end

    # discrete callback to update Es to ensure continuity of η
    grade_change_condition(s, x, integrator) = x in integrator.p.eetcprob.track.x_gradient

    function grade_change_aff!(integrator)
        params = integrator.p
        track = params.eetcprob.track
        train = params.eetcprob.train
        Es = params.Es
        η₀ = calculate_η(integrator.u, integrator.tprev, params)   # current η
        g_old = g(track, integrator.t - 0.1)
        g_new = g(track, integrator.t + 0.1)
        if params.current_phase != MaxB
            E_old = last(Es)
            E_new = (g_new - g_old) * η₀ + E_old
            push!(Es, E_new)
        else
            F_old = last(Es)
            F_new = (g_new - g_old) * (η₀ - train.ρ + 1.) + F_old
            push!(Es, F_new)
        end
    end

    odeprob = ODEProblem(_odefun, s0, xspan, simparams)

    # setup callbacks
    vector_cont_cb = VectorContinuousCallback(
        mode_switch_condition,
        mode_switch_affect_pos!,
        mode_switch_affect_neg!,
        VECTOR_CB_LENGTH,
        save_positions = (false, false),
        rootfind = SciMLBase.LeftRootFind
    )

    # saving callback
    saved_vals = SavedValues(Float64, Vector{Union{Float64, Mode}})  # time type, η type
    function saving_func(s, x, integrator)
        if x in integrator.p.eetcprob.track.x_gradient
            [calculate_η(s, integrator.tprev, integrator.p), integrator.p.current_phase]
        else
            [calculate_η(s, x, integrator.p), integrator.p.current_phase]
        end
    end
    saving_cb = SavingCallback(saving_func, saved_vals)

    # discrete callback
    grade_cb = DiscreteCallback(grade_change_condition, grade_change_aff!,
        save_positions = (false, false))

    cb_set = CallbackSet(vector_cont_cb, saving_cb, grade_cb)

    odesol = OrdinaryDiffEq.solve(odeprob, AutoVern7(Rodas5()),
        callback = cb_set,
        dtmax = 3.,
        tstops = track.x_gradient,
        d_discontinuities = track.x_gradient,)

    η_and_mode = saved_vals.saveval |> stack
    η = η_and_mode[1,:][:]
    η = Vector{Float64}(η)
    modes = η_and_mode[2,:][:]
    modes = Vector{Mode}(modes)

    x_mode_switch = [0.; odesol.t[findall(k -> modes[k] != modes[k-1], 2:length(modes)) .+ 1]]
    modes_sequence = modes[[1; findall(k -> modes[k] != modes[k-1], 2:length(modes)) .+ 1]]

    control = function (x)
        current_phase = modes_sequence[searchsortedlast(x_mode_switch, x)]
        if current_phase == MaxP
            simparams.eetcprob.train.U̅(odesol(x)[2])
        elseif current_phase == Coast
            zero(x)
        elseif current_phase == MaxB
            simparams.eetcprob.train.U̲(odesol(x)[2])
        end
    end

    otc_sol = OTCSolution(
        odesol,
        x_mode_switch,
        modes_sequence,
        control,
        η    
    )
end

# simulate
otc_sol = simulate_regular_forward(simparams, xspan, s0)
odesol = otc_sol.odesol
η = otc_sol.η

# plot(odesol.t, modes)
plot(odesol.t, η)
# plot(η, odesol[2,:])

## Profiling
function prof_g(n, simparams, xspan, s0)
    for _ in 1:n
        simulate_regular_forward(simparams, xspan, s0)
    end
end

@profview_allocs prof_g(100, simparams, xspan, s0)