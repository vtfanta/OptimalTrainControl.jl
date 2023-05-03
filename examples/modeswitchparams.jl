import ForwardDiff: derivative
using DifferentialEquations
using NumericalIntegration
using Plots
using OptimalTrainControl
using Roots

@enum ControlModes begin
    MaxP
    HoldP
    Coast
    MaxB
end

"""
Return solution for given speed V and given HoldP length.
"""
function __solveflat(v₀, vf, V, track::Track, res::Resistance, ρ, 
    maxcontrol = x -> 1. / max(5., x), mincontrol = x -> -1. / max(5., x), holdPstart = nothing, holdPlength = nothing)
    function mycontrol(u, p, x)
        v = u[2]
        if p == MaxP
            maxcontrol(v)
        elseif p == HoldP
            resistance(res, v) - getgradientacceleration(track, x)
        elseif p == Coast
            0.0
        elseif p == MaxB
            mincontrol(v)
        end
    end
    
    function odefun!(du, u, p, x)
        t, v, μ₂ = u
    
        η = μ₂ / v - 1.
        ζ = η + 1 - ρ
    
        π₁ = η > 0. ? η : 0.
        π₂ = ζ ≥ 0. ? 0. : -ζ

        # @show v, μ₂, x, mycontrol(v, μ₂, x)
    
        du[1] = inv(v)
        du[2] = (mycontrol(u, p, x) - resistance(res, v) + 
            getgradientacceleration(track, x)) * inv(v)
        du[3] = -ψ(res, V) / v^2 + μ₂ * du[2] / v + 
            μ₂ * derivative(x -> resistance(res, x), v) / v - 
            π₁ * derivative(x -> maxcontrol(x), v) 
            π₂ * derivative(x -> mincontrol(x), v)
    end

    # Calculate the initial value of μ₂
    # Assuming that v₀ < V!!
    μ₂0 = v₀ * ( (E(res, V, v₀) - E(res, V, V)) / 
        (maxcontrol(v₀) - resistance(res, v₀) + getgradientacceleration(track, start(track))) + 1.)

    function condition(out, u, t, _)
        # Terminating at reaching the final speed
        out[1] = u[2] - vf
        # Reaching cruising speed, switch to HoldP or Coast
        out[2] = u[2] - V
        # Reaching end of the HoldP phase
        out[3] = isnothing(holdPstart) ? 1 : t - (holdPstart + holdPlength)
        # Start the braking phase
        out[4] = ρ * u[2] - u[3]
    end
    function affect!(int, idx)
        # @show int.t
        if idx == 1
            terminate!(int)
        elseif !isnothing(holdPstart) && idx == 2
            @debug "Hitting V, to singular HoldP"
            int.p = HoldP
            holdPstart = int.t
        elseif isnothing(holdPstart) && idx == 2
            @debug "Hitting V, to regular Coast"
            # @show int.u[2], int.u[3]
            int.p = Coast
        elseif idx == 3
            @debug "Ending HoldP, to regular Coast"
            int.p = Coast
        elseif idx == 4
            @debug "Starting MaxB"
            int.p = MaxB
        end
    end
    cb = VectorContinuousCallback(condition, affect!, 4)

    prob = ODEProblem(odefun!, [0., v₀, μ₂0], (start(track), Inf), MaxP)

    if isnothing(holdPstart)
        sol = solve(prob, alg_hints = [:stiff], callback = cb)
    else
        sol = solve(prob, alg_hints = [:stiff], callback = cb, 
            tstops = [holdPstart, holdPstart + holdPlength])
    end

    x = sol.t
    t = sol[1,:]
    v = sol[2,:]
    μ₂ = sol[3,:]

    η = μ₂ ./ v .- 1.

    x, t, v, η
end

"Return solution for given speed V, HoldP phase length is adjusted automatically."
function _solveflat(v₀, vf, track, res, ρ, V, 
    maxcontrol = x -> 1. / max(5., x), mincontrol = x -> -1. / max(5., x))
    
    x, t, v, η = __solveflat(v₀, vf, V, track, res, ρ, maxcontrol, mincontrol)

    holdPstart = x[argmax(v)]
    holdPlength = finish(track) - x[end]
    if holdPlength ≤ 0
        holdPstart, holdPlength = nothing, nothing
    end

    # display(plot(x,v))
    @show x[end]
    @info "Adjusting length of the cruising phase to $holdPlength m."

    x, t, v, η = __solveflat(v₀, vf, V, track, res, ρ, maxcontrol, mincontrol, holdPstart, holdPlength)
end

"Return solution for given total time T."
function solveflat(v₀, vf, track, res, ρ, T, 
    maxcontrol = x -> 1. / max(5., x), mincontrol = x -> -1. / max(5., x))

    mintimemodel = AlbrechtModel(res, maxcontrol, mincontrol, 1e3)
    mintimescenario = MinimalTimeScenario(mintimemodel, track, 9.81, [v₀], [vf])
    calculatecontrol!(mintimescenario)
    sol = play(mintimescenario)
    maxV = maximum(sol[1,:]) - 1
    minT = integrate(sol.t, inv.(sol[1,:]))

    if T ≤ minT
        error("The requested time $T s is not feasible. Best achievable is $minT s.")
    end

    @show minT

    if abs(T - minT) / minT ≤ 0.1
        @warn "Requested time is close to the minimal one. Results may be wrong."
    end

    function timedifference(candidateV)
        x, t, _, _ = _solveflat(v₀, vf, track, res, ρ, candidateV, maxcontrol, mincontrol)
        if x[end] ≥ finish(track) + 10.0
            return x[end] - finish(track)
        end
        t[end] - T
    end

    Vs = collect(finish(track)/(2T):2:maxV)

    display(scatter(Vs, timedifference.(Vs)))

    # timedifference(finish(track) / T)
    # V = find_zero(v -> timedifference(v), finish(track)/T)
    @info "Travel speed found: $V m/s."

    x, t, v, η = _solveflat(v₀, vf, track, res, ρ, V, maxcontrol, mincontrol)
end

v₀, vf = 1.0, 1.0

ρ = 0.5

V = 19.8

T = 1050

mytrack = FlatTrack(20e3)

myres = DavisResistance(1e-2, 0., 1.5e-5)

x, t, v, η = solveflat(v₀, vf, mytrack, myres, ρ, T)

# display(plot(η, v, xlabel = "η", ylabel = "v (m/s)", legend = false))
plot(x, v, xlabel = "Distance (m)", ylabel = "v (m/s)", legend = false)