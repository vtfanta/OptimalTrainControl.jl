using OrdinaryDiffEq
using Roots

module LinkPortConditionModule
    struct LinkPortCondition{T<:Real}
        end_port_start::T
        end_port_finish::T
        end_port_speed::T
        end_port_η::T
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

# simulate with regular modes with terminate! callback when target speed on the target link is reached
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

function root_f_start(E1, eetcprob, V, W, cond::LinkPortConditionModule.LinkPortCondition{T}) where {T<:Real}
    Es = [E1]
    
    if eetcprob.initial_speed > V
        initial_mode = Coast
    else    # eetcprob.initial_speed ≤ V
        initial_mode = MaxP
    end
    
    simparams = EETCSimParams(eetcprob, V, W, Es, initial_mode)
    otc_sol = simulate_link_forward(zero(T), cond, simparams, SA[zero(T), eetcprob.initial_speed]);
    
    η = otc_sol.η
    retcode = otc_sol.odesol.retcode

    if retcode == SciMLBase.ReturnCode.Terminated && otc_sol.odesol[2,end] ≤ 0.5    # terminated, low speed
        return +Inf        
    elseif retcode == SciMLBase.ReturnCode.Terminated && abs(otc_sol.odesol[2,end] - V) ≤ 1e-2    # terminated, reached target speed
        return η[end] - cond.end_port_η
    elseif retcode == SciMLBase.ReturnCode.Success  # reached end of target port
        return +Inf
    end
end

# link starting port means finding such Es[1] that η has proper value at reaching target port with right speed
# this is done via root finding
function link_start(port2::Port{T}, simparams::EETCSimParams{T,F1,F2}; atol=1e-6) where {T<:Real,F1,F2}
    target_speed = port2.speed

    if port2.mode == HoldP || port2.mode == HoldP_SL
        target_η = zero(T)
    elseif port2.mode == HoldR || port2.mode == HoldR_SL
        target_η = simparams.eetcprob.train.ρ - one(T)
    end

    port_condition = LinkPortConditionModule.LinkPortCondition(
        port2.start,
        port2.finish,
        target_speed,
        target_η
    )
    E1 = Roots.find_zero(
        e -> root_f_start(e, simparams.eetcprob, simparams.V, simparams.W, port_condition),
        [-30., 30.], A42(),
        atol = atol
    )
    simulate_link_forward(
        zero(T),
        port_condition, 
        EETCSimParams(simparams.eetcprob, simparams.V, simparams.W, [E1], simparams.current_phase),
        SA[zero(T), simparams.eetcprob.initial_speed])
end

# choose proper function to link two ports
function link(port1::Port{T}, port2::Port{T}, simparams::EETCSimParams{T,F1,F2}) where {T,F1,F2}
    if isinf(port1.start)   # starting port has starting point at -∞
        # link starting port
        return link_start(port2, simparams)
    elseif isinf(port2.finish)   # ending port has finishing point at +∞
        # TODO
    else    # link two interior ports
        # TODO
    end
end