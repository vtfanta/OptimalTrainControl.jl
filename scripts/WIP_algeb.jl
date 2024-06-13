using DiffEqCallbacks
using OptimalTrainControl
using OrdinaryDiffEq
using Plots
using Roots
using StaticArrays

function make_cond_power2coast(V)
    return function cond_power2coast(s, x, int)
        η = calculate_η(s, x, int, V)
        η
    end
end

function make_cond_coast2brake(V)
    return function cond_coast2brake(s, x, int)
        η = calculate_η(s, x, int, V)
        η - (int.p.train.ρ - 1.)
    end
end

aff_2coast!(int) = int.p.current_phase = Coast

aff_2brake!(int) = int.p.current_phase == Coast ? int.p.current_phase = MaxB : nothing

aff_2power!(int) = int.p.current_phase == Coast ? int.p.current_phase = MaxP : nothing

function make_update_Es(first_η)
    last_η = first_η
    return function update_Es!(s, x, int)
        # @show x, last_η
        if x in int.p.track.x_gradient
            g_prev = g(int.p.track, x - 1e-2)
            g_next = g(int.p.track, x + 1e-2)
            if int.p.current_phase == MaxP || int.p.current_phase == Coast
                push!(int.p.Es, (g_next - g_prev) * last_η + last(int.p.Es))
            elseif int.p.current_phase == MaxB
                # this is actually updating F_i for ζ
                push!(int.p.Es, (g_next - g_prev) * (last_η + 1. - int.p.train.ρ)) + last(int.p.Es)
            end
        else
            last_η = calculate_η(s, x, int, V)
        end
        # @show last_η
    end
end

function _odefun(s::A, p::EETCProblem, x::T) where {T<:Real, A<:AbstractArray{T,1}}
    t, v = s

    if p.current_phase == MaxP
        u = p.train.U̅(v)
    elseif p.current_phase == Coast
        u = 0.
    elseif p.current_phase == MaxB
        u = p.train.U̲(v)
    end
    # @show u

    ds1 = 1/v
    ds2 = (u - OptimalTrainControl.r(p.train, v) + OptimalTrainControl.g(p.track, x)) / v
    SA[ds1, ds2]
end

prob = EETCProblem(;
    train,
    track = track2,
    current_phase = MaxP,
    initial_speed = V,
    T = 1e3,
    Es = [-E(V, V)]
)

xspan = (1942.77, 3100.)
s0 = SA[0., V]
odeprob = ODEProblem(_odefun, s0, xspan, prob;
    tstops = track2.x_gradient)

# low speed termination (0.5 m/s)
lowspeed_cond(s, x, int) = s[2] - 0.5
lowspeed_aff(int) = terminate!(int)
lowspeed_cb = ContinuousCallback(lowspeed_cond, lowspeed_aff)

# end after end of the steep section when hitting V
targetspeed_cond(s, x, int) = s[2] - V
targetspeed_aff(int) = terminate!(int)
targetspeed_cb = ContinuousCallback(targetspeed_cond, targetspeed_aff; affect_neg! = nothing)

# saving η (ζ)
function calculate_η(s, x, int, V)
    v = s[2]
    
    ψ(v) = (int.p.train.r[2] + 2int.p.train.r[3]*v) * v^2
    W = Roots.find_zero(v -> -ψ(V) + int.p.train.ρ*ψ(v), V)
    E(V, v) = ψ(V)/v + OptimalTrainControl.r(int.p.train, v)

    if int.p.current_phase == MaxP
        val = (E(V, v) + last(int.p.Es)) / 
            (int.p.train.U̅(v) - OptimalTrainControl.r(int.p.train, v) + g(int.p.track, x))
        # @show val
    elseif int.p.current_phase == Coast
        val = (E(V, v) + last(int.p.Es)) / 
            (- OptimalTrainControl.r(int.p.train, v) + g(int.p.track, x))
    elseif int.p.current_phase == MaxB
        # have to calculate ζ, the Fs are also in the EETC.Es array
        val = (int.p.train.ρ*E(W,v) + last(int.p.Es)) /
            (int.p.train.U̲(v) - OptimalTrainControl.r(int.p.train, v) + g(int.p.train, x))
        # η = ζ + ρ - 1
        val += int.p.train.ρ - 1.
    end
    return val
end
saved_vals = SavedValues(Float64, Float64)  # time type, savedval type
saving_cb = SavingCallback((s,x,int) -> calculate_η(s,x,int,V), saved_vals)

powercoast_cb = ContinuousCallback(make_cond_power2coast(V),aff_2power!,
    affect_neg! = aff_2coast!)

coastbrake_cb = ContinuousCallback(make_cond_coast2brake(V), aff_2coast!,
    affect_neg! = aff_2brake!)

tstops = copy(track2.x_gradient)
push!(tstops, (clamp.(track2.x_gradient .- 1.0, 0., Inf))...) 

updateEs_cb = FunctionCallingCallback(make_update_Es(0.); funcat = tstops[tstops .> 0.], func_start = false)

callbacks = CallbackSet(lowspeed_cb, targetspeed_cb, updateEs_cb, saving_cb, coastbrake_cb)

odesol = OrdinaryDiffEq.solve(odeprob, Tsit5();
    callback = callbacks,
    d_discontinuities = track2.x_gradient,
    tstops,
    dtmax = 100.)
plot(odesol.t, odesol[2,:])
# @show prob.Es
η = saved_vals.saveval
push!(η, 0.)
plot(η)

# Es = [-E(V, V)]
# for k in eachindex(odesol.t)
#     v = odesol[2,k]
#     x = odesol.t[k]

#     if x in track2.x_gradient
#         # @show x
#         # η[k] = (E(V, v) + last(Es)) / (train.U̅(v) - OptimalTrainControl.r(train, v) + g(track2, x-0.1))
#         # @show η[k]
#         push!(Es, (g(track2, x+1.0) - g(track2, x-1.0)) * η[k-1] + last(Es))
#         # @show (E(V, v) + last(Es)) / (train.U̅(v) - OptimalTrainControl.r(train, v) + g(track2, x))
#         η[k] = (E(V, v) + last(Es)) / (train.U̅(v) - OptimalTrainControl.r(train, v) + g(track2, x))
#     else
#         η[k] = (E(V, v) + last(Es)) / (train.U̅(v) - OptimalTrainControl.r(train, v) + g(track2, x))
#     end
# end
# @show η[end]
# plot(odesol.t, η)

# plot(η, odesol[2,:])
# xlims!(-5e-3, 5e-3)
# ylims!(0, 12)

##

# What I need:
#   - callback for transition from MaxP to Coast and vice versa with _neg
#   - callback for transition from Coast to MaxB and vice versa with _neg
#   x callback to automatically push! next E at the points of grade change
#   x callback for saving adjoint variables at every timestep