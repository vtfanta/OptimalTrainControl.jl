using DiffEqCallbacks
using OptimalTrainControl
using OrdinaryDiffEq
using Plots
using Roots
using StaticArrays

function make_cond_power2coast(p::EETCProblem, V)
    return function cond_power2coast(s, x, int)
        η = calculate_η(s, x, int, V)
        η
    end
end

function make_cond_coast2brake(p::EETCProblem, V)
    return function cond_coast2brake(s, x, int)
        η = calculate_η(s, x, int, V)
        η - (int.p.train.ρ - 1.)
    end
end

aff_2coast!(int) = int.p.current_phase = Coast

function make_aff_2brake!(p::EETCProblem, V)

    function aff_2brake!(int)
        int.p.current_phase == Coast ? int.p.current_phase = MaxB : nothing

        # add the last F to the E array
        # ζ = 0 -> F = -ρ * E(W, v)
        push!(int.p.Es, )
    end
end

aff_2power!(int) = int.p.current_phase == Coast ? int.p.current_phase = MaxP : nothing

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

# saving η (ζ shifted)
function calculate_η(s, x, int, V)
    v = s[2]
    
    if int.p.current_phase == MaxP
        val = (E(int.p.train, V, v) + last(int.p.Es)) / 
            (int.p.train.U̅(v) - OptimalTrainControl.r(int.p.train, v) + g(int.p.track, x))
        # @show val
    elseif int.p.current_phase == Coast
        val = (E(int.p.train, V, v) + last(int.p.Es)) / 
            (- OptimalTrainControl.r(int.p.train, v) + g(int.p.track, x))
    elseif int.p.current_phase == MaxB
        # have to calculate ζ, the Fs are also in the EETC.Es array
        val = (int.p.train.ρ*E(int.p.train, W,v) + last(int.p.Es)) /
            (int.p.train.U̲(v) - OptimalTrainControl.r(int.p.train, v) + g(int.p.track, x))
        # η = ζ + ρ - 1
        val += int.p.train.ρ - 1.
    end
    return val
end

# calculate E coefficients for individual gradient transitions
function make_update_Es(V)
    return function update_Es!(s, x, int)
        p = int.p
        η_current = calculate_η(s, x-0.01, int, V) # just before gradient change
        g_prev = g(p.track, x-0.01)
        g_next = g(p.track, x+0.01)
        if p.current_phase == MaxP || p.current_phase == Coast  # save E
            push!(p.Es, (g_next - g_prev) * η_current + last(p.Es))
        elseif p.current_phase == MaxB  # save F
            push!(p.Es, (g_next - g_prev) * (η_current + 1 - p.train.ρ) + last(p.Es))
        end
    end
end

train = Train(
    v -> 1/v,
    v -> -1/v,
    (1e-2, 0., 1.5e-5),
    0.6
)

uphill_track = Track(;
    length = 5e3,
    x_gradient = [0., 2e3, 3e3],
    gradient = [0., 35/1e3, 0.]
)

V = 10.

prob = EETCProblem(;
    train,
    track = uphill_track,
    current_phase = MaxP,
    initial_speed = V,
    T = 1e3,
    Es = [-E(train, V, V)]  # To ensure η = 0 at the beginning (HoldP)
)

s0 = SA[0., V]
xspan = (1910, 2110)
odeprob = ODEProblem(_odefun, s0, xspan, prob;
    tstops = uphill_track.x_gradient)

# low speed termination (0.5 m/s)
lowspeed_cond(s, x, int) = s[2] - 0.5
lowspeed_aff(int) = terminate!(int)
lowspeed_cb = ContinuousCallback(lowspeed_cond, lowspeed_aff)

# end after end of the steep section when hitting V
targetspeed_cond(s, x, int) = s[2] - V
targetspeed_aff(int) = terminate!(int)
targetspeed_cb = ContinuousCallback(targetspeed_cond, targetspeed_aff; affect_neg! = nothing)

# saving η
saved_vals = SavedValues(Float64, Float64)  # time type, savedval type
saving_cb = SavingCallback((s,x,int) -> calculate_η(s,x,int,V), saved_vals)

powercoast_cb = ContinuousCallback(make_cond_power2coast(V),aff_2power!,
    affect_neg! = aff_2coast!)

coastbrake_cb = ContinuousCallback(make_cond_coast2brake(V), aff_2coast!,
    affect_neg! = aff_2brake!)

tstops = copy(uphill_track.x_gradient)
# push!(tstops, (clamp.(track2.x_gradient .- 1.0, 0., Inf))...) 

updateEs_cb = FunctionCallingCallback(make_update_Es(0., V); funcat = uphill_track.x_gradient[2:end], func_start = false)

callbacks = CallbackSet(lowspeed_cb, targetspeed_cb, updateEs_cb, saving_cb)

odesol = OrdinaryDiffEq.solve(odeprob, Tsit5();
    callback = callbacks,
    d_discontinuities = uphill_track.x_gradient,
    tstops,
    dtmax = 10.)
plot(odesol.t, odesol[2,:])