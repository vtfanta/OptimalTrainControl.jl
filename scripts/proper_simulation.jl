using DiffEqCallbacks
using OptimalTrainControl
using OrdinaryDiffEq
using Plots
using StaticArrays

# use this script to finetune the process of running simulation of regular phases (MaxP, Coast, MaxB)
# - ensure proper transition between phases
# - ensure proper calculation and saving of η
# - ensure proper calculations of Es (necessary for η)

module ParamsWrapper
    using OptimalTrainControl
    mutable struct EETCSimParams{T<:Real}
        eetcprob::EETCProblem
        V::T
        W::T
        Es::Vector{T}
        current_phase::Mode
    end
end

EETCSimParams = ParamsWrapper.EETCSimParams # TODO remove this wrapper, WARNING

# define the track
track = Track(
    length = 5e3,
    altitude = 10.,
    x_gradient = [0., 2e3, 3e3, 4e3],
    gradient = [0., 35/1e3, 0., -35/1e3],
)

# define the train
train = Train(
    v -> 1/v,
    v -> -1/v,
    (1e-2, 0., 1.5e-5),
    0.6
)

# guess initial Es[1]; during solution of EETC problem, this needs to be optimized for
Es1 = -0.0185

# define the EETC problem
total_time = 500.
prob = EETCProblem(total_time, train, track)
prob.Es = [Es1]             # TODO remove from EETCProblem struct
prob.current_phase = MaxP   # TODO remove from EETCProblem struct
prob.initial_speed = 1.

simparams = EETCSimParams(prob, 10., OptimalTrainControl.calculate_W(prob, 10.), [Es1], Coast)


"""
In-place version of the right-hand side of the ODE system for the EETC problem.

Params are in a EETCSimParams struct.
"""
function _odefun(s::A, p::EETCSimParams{T}, x::T) where {T<:Real, A<:AbstractArray{T,1}}
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
    if phase == MaxP
        val = (E(train, V, v) + last(Es)) / 
            (train.U̅(v) - OptimalTrainControl.r(train, v) + g(track, x))
    elseif phase == Coast
        val = (E(train, V, v) + last(Es)) / 
            (- OptimalTrainControl.r(train, v) + g(track, x))
    elseif phase == MaxB
        # have to calculate ζ, the Fs are also in the EETC.Es array
        # TODO can E be found such that I dont have to store F instead?
        val = (train.ρ*E(train, W, v) + last(Es)) /
            (train.U̲(v) - OptimalTrainControl.r(train, v) + g(track, x))
        # η = ζ + ρ - 1
        val += train.ρ - 1.
    end
    return val
end

function mode_switch_condition(out, s, x, integrator)
    _, v = s
    p = integrator.p
    η = calculate_η(s, x, p)
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

# setup simulation
s0 = SA[0., 10.]
xspan = (1980., 2300)
odeprob = ODEProblem(_odefun, s0, xspan, simparams)

# setup callbacks
const VECTOR_CB_LENGTH = 3
vector_cont_cb = VectorContinuousCallback(
    mode_switch_condition,
    mode_switch_affect_pos!,
    mode_switch_affect_neg!,
    VECTOR_CB_LENGTH,
    save_positions = (false, false),
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
modes = η_and_mode[2,:][:]

plot(odesol.t, modes .|> Int)
plot(odesol.t, η)
plot(η, odesol[2,:])