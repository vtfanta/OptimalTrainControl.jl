using OrdinaryDiffEq

module LinkPortConditionModule
    struct LinkPortCondition{T<:Real}
        end_port_start::T
        end_port_finish::T
        end_port_speed::T
    end
end
function (c::LinkPortConditionModule.LinkPortCondition)(s, x::T, integrator) where {T<:Real}
    s[2] - c.end_port_speed
end
function (c::LinkPortConditionModule.LinkPortCondition)(integrator)
    if c.end_port_start ≤ integrator.t ≤ c.end_port_finish
        terminate!(integrator)
    end
end

function simulate_link_forward(xstart::T, cond::LinkPortConditionModule.LinkPortCondition{T}, sparams::EETCSimParams{T,F1,F2}, s0::SVector{2, T}) where {T<:Real, F1, F2}
    xspan = (xstart, cond.end_port_finish)
    
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
    link_cont_cb = ContinuousCallback(cond, cond, save_positions = (false, false))

    # saving callback
    saved_vals = DiffEqCallbacks.SavedValues(eltype(xspan), Vector{Union{eltype(xspan), Mode}})  # time type, η type
    saving_cb = SavingCallback(saving_func, saved_vals)

    # discrete callback
    grade_cb = DiscreteCallback(grade_change_condition, grade_change_aff!,
        save_positions = (false, false))

    cb_set = CallbackSet(vector_cont_cb, saving_cb, grade_cb, link_cont_cb)

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