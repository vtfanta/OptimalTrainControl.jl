using OptimalTrainControl
using OrdinaryDiffEq
using DiffEqCallbacks

# use this script to finetune the process of running simulation of regular phases (MaxP, Coast, MaxB)
# - ensure proper transition between phases
# - ensure proper calculation and saving of η
# - ensure proper calculations of Es (necessary for η)

const EETCPROB_INDEX = 1
const 

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
Es1 = -0.0285

# define the EETC problem
total_time = 500.
prob = EETCProblem(total_time, train, track)
prob.Es = [Es1]
prob.current_phase = MaxP
prob.initial_speed = 1.

"""
In-place version of the right-hand side of the ODE system for the EETC problem.

Params are in named tuple (EETCProblem, V, W, Es, current_phase)
"""
function _odefun!(ds::A, s::A, p::@NamedTuple{eetcprob::EETCProblem, V::Real, W::Real, Es::Vector{T}, current_phase::Mode}, x::T) where {T<:Real, A<:AbstractArray{T,1}}
    _, v = s

    eetcprob, _, _, _, phase = p
    if phase == MaxP
        u = eetcprob.train.U̅(v)
    elseif phase == Coast
        u = 0.
    elseif phase == MaxB
        u = eetcprob.train.U̲(v)
    end
    # @show u

    ds[1] = 1/v
    ds[2] = (u - OptimalTrainControl.r(eetcprob.train, v) + OptimalTrainControl.g(eetcprob.track, x)) / v
end 

function calculate_η(s, x, p::@NamedTuple{eetcprob::EETCProblem, V::Real, W::Real, Es::Vector{T}, current_phase::Mode})
    v = s[2]
    
    eetcprob, V, W, Es, phase = p
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
    out[2] = η - (p[1].train.ρ - 1.)    # boundary between Coast and MaxB
    # TODO should I add terminating condition?
    out[3] = v - 0.1    # terminate when speed is too low
end

function mode_switch_affect_pos!(integrator, index) # gets applied when crossing 0 from negative to positive
    if index == 1   # Coast to MaxP
        integrator.p[end] = MaxP
    elseif index == 2   # MaxB to Coast
        integrator.p[end] = Coast
        # append to Es since switching from ζ to η
        η₀ = calculate_η(integrator.u, integrator.t, integrator.p)
        v = integrator.u[2]
        Es = integrator.p
        eetcprob = integrator.p.eetcprob
        train = eetcprob.train
        track = eetcprob.track
        
        newE = η₀ * (-r(train, v) + g(track, integrator.t)) - E(train, integrator.p.V, v)
        push!(Es, newE)
    elseif index == 3   # terminate because of low speed
        terminate!(integrator)
    end
end

function mode_switch_affect_neg!(integrator, index) # gets applied when crossing 0 from positive to negative
    if index == 1   # MaxP to Coast
        integrator.p[end] = Coast
    elseif index == 2   # Coast to MaxB 
        η₀ = calculate_η(integrator.u, integrator.t, integrator.p)
        v = integrator.u[2]
        Es = integrator.p
        eetcprob = integrator.p.eetcprob
        train = eetcprob.train
        track = eetcprob.track

        integrator.p[end] = MaxB

        # append to Es since switching from η to ζ (it's actually F)
        newF = η₀ * (train.U̲(v) - r(train, v) + g(track, integrator.t)) - train.ρ * E(train, integrator.p.W, v)
        push!(Es, newF)
    end
end