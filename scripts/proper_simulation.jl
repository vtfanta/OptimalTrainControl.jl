using OptimalTrainControl
using OrdinaryDiffEq
using DiffEqCallbacks

# use this script to fine tune the process of running simulation of regular phases (MaxP, Coast, MaxB)
# - ensure proper transition between phases
# - ensure proper calculation and saving of η
# - ensure proper calculations of Es (necessary for η)

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

Params are in tuple (EETCProblem, V, W, Es, current_phase)
"""
function _odefun!(ds::A, s::A, p::Tuple{EETCProblem, Real, Real, Vector{T}, Mode}, x::T) where {T<:Real, A<:AbstractArray{T,1}}
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

function calculate_η(s, x, p::Tuple{EETCProblem, Real, Real, Vector{T}, Mode})
    v = s[2]
    
    eetcprob, V, W, Es, phase = p
    train = eetcprob.train
    track = eetcprob.track
    if phase == MaxP
        val = (E(train, V, v) + last(Es)) / 
            (train.U̅(v) - OptimalTrainControl.r(train, v) + g(track, x))
        # @show val, x, p.current_phase
    elseif phase == Coast
        val = (E(train, V, v) + last(Es)) / 
            (- OptimalTrainControl.r(train, v) + g(track, x))
    elseif phase == MaxB
        # have to calculate ζ, the Fs are also in the EETC.Es array
        # TODO finish this, can E be found such that I dont have to store F instead?
        val = (p.train.ρ*E(p.train, W,v) + last(p.Es)) /
            (p.train.U̲(v) - OptimalTrainControl.r(p.train, v) + g(p.track, x))
        # η = ζ + ρ - 1
        val += p.train.ρ - 1.
    end
    return val
end

function mode_switch_condition(out, s, x, integrator)
    
end
