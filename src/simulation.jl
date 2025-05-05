using DiffEqCallbacks

export calculate_η, simulate_regular_forward

const VECTOR_CB_LENGTH = 3

"""
Out-of-place version of the right-hand side of the ODE system for the EETC problem.

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
        val = (train.ρ*E(train, W, v) + Es[end]) /
            (train.U̲(v) - OptimalTrainControl.r(train, v) + g(track, x))
        # η = ζ + ρ - 1
        val += train.ρ - 1.
    end
    return val
end

function mode_switch_condition(out::A, s::B, x::T, integrator) where {T<:Real,A<:AbstractArray{T,1}, B<:AbstractArray{T,1}}
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

function saving_func(s, x, integrator)
    if x in integrator.p.eetcprob.track.x_gradient
        [calculate_η(s, integrator.tprev, integrator.p), integrator.p.current_phase]
    else
        [calculate_η(s, x, integrator.p), integrator.p.current_phase]
    end
end

function simulate_regular_forward(sparams::EETCSimParams{T,F1,F2}, xspan::Tuple{T, T}, s0::SVector{2, T}) where {T<:Real, F1, F2}
    
    # to prevent mutating
    simparams = deepcopy(sparams)

    odeprob = ODEProblem{false}(_odefun, s0, xspan, simparams)

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
    saved_vals = DiffEqCallbacks.SavedValues(eltype(xspan), Vector{Union{eltype(xspan), Mode}})  # time type, η type
    saving_cb = SavingCallback(saving_func, saved_vals)

    # discrete callback
    grade_cb = DiscreteCallback(grade_change_condition, grade_change_aff!,
        save_positions = (false, false))

    cb_set = CallbackSet(vector_cont_cb, saving_cb, grade_cb)

    odesol = OrdinaryDiffEq.solve(odeprob, AutoVern7(Rodas5()),
        callback = cb_set,
        dtmax = 3.,
        tstops = simparams.eetcprob.track.x_gradient,
        d_discontinuities = simparams.eetcprob.track.x_gradient)

    η_and_mode = saved_vals.saveval

    sol_length = length(odesol.t)
    η = Vector{eltype(xspan)}(undef, sol_length)
    modes = Vector{Mode}(undef, sol_length)
    for i in 1:sol_length
        η[i] = η_and_mode[i][1]
        modes[i] = η_and_mode[i][2]
    end

    indices = [1; findall(k -> modes[k] != modes[k-1], 2:length(modes)) .+ 1]
    x_mode_switch = odesol.t[indices]
    x_mode_switch = Vector{eltype(xspan)}(x_mode_switch)
    modes_sequence = modes[indices]

    control = OptimalTrainControl.create_concrete_control_function(
        x_mode_switch,
        modes_sequence,
        odesol,
        simparams.eetcprob.train
    )

    OTCSolution(
        odesol,
        x_mode_switch,
        modes_sequence,
        control,
        η    
    )
end