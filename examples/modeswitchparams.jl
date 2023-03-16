import ForwardDiff: derivative
using DifferentialEquations
using NumericalIntegration
using Plots
using RailDynamics
using Roots

@enum ControlModes begin
    MaxP
    HoldP
    Coast
    MaxB
end

function __solveflat(v₀, vf, V, track::Track, res::Resistance, ρ, justswitchingpoints = false,
    maxcontrol = x -> 1. / max(5., x), mincontrol = x -> -1. / max(5., x), holdPlength = 0.0)
    holdPstart = nothing
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
        out[3] = isnothing(holdPstart) ? 1 : t - (holdPlength + holdPstart)
        # Start the braking phase
        out[4] = ρ * u[2] - u[3]
    end
    function affect!(int, idx)
        @show int.t
        if idx == 1
            terminate!(int)
        elseif idx == 2 && holdPlength > 0.1
            @info "Hitting V, to singular HoldP"
            int.p = HoldP
            holdPstart = int.t
        elseif idx == 2 && holdPlength ≈ 0.0
            @info "Hitting V, to regular Coast"
            # @show int.u[2], int.u[3]
            int.p = Coast
        elseif idx == 3
            @info "Ending HoldP, to regular Coast"
            int.p = Coast
        elseif idx == 4
            @info "Starting MaxB"
            int.p = MaxB
        end
    end
    cb = VectorContinuousCallback(condition, affect!, 4)

    prob = ODEProblem(odefun!, [0., v₀, μ₂0], (start(track), Inf), MaxP)


    sol = solve(prob, alg_hints = [:stiff], callback = cb)
    x = sol.t
    t = sol[1,:]
    v = sol[2,:]
    μ₂ = sol[3,:]

    η = μ₂ ./ v .- 1.

    if justswitchingpoints
        @warn "Not implemented yet/"
    else
        x, t, v, η
    end
end

function _solveflat(v₀, vf, track, res, ρ, V, justswitchingpoints = false, 
    maxcontrol = x -> 1. / max(5., x), mincontrol = x -> -1. / max(5., x))
    
    x, t, v, η = __solveflat(v₀, vf, V, track, res, ρ, justswitchingpoints, maxcontrol, mincontrol)


    holdPlength = finish(track) - x[end]
    # @show holdPlength, x[end]
    
    cnt = 0
    while true
        x, t, v, η = __solveflat(v₀, vf, V, track, res, ρ, justswitchingpoints, maxcontrol, mincontrol, holdPlength)
        Δend = x[end] - finish(track)
        @show finish(track), Δend
        if abs(Δend) < 10.0 || cnt >= 2
            # display(plot(x,v))
            return x, t, v, η
        else
            @info "Adjusting the HoldP length"
            holdPlength -= 0.8 * Δend
            display(plot(x,v))
            cnt += 1
        end
    end
end

function solveflat(v₀, vf, track, res, ρ, T, justswitchingpoints = false, 
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


    function timedifference(candidateV)
        _, t, _, _ = _solveflat(v₀, vf, track, res, ρ, candidateV, false, maxcontrol, mincontrol)
        t[end] - T
    end

    Vs = collect(max(v₀, vf):2:maxV)

    display(scatter(Vs, timedifference.(Vs)))

    timedifference(finish(track) / T)
end

v₀, vf = 1.0, 1.0

ρ = 0.5

V = 20

T = 1000

mytrack = FlatTrack(25e3)

myres = DavisResistance(1e-2, 0., 1.5e-5)

x, t, v, η = _solveflat(v₀, vf, mytrack, myres, ρ, V)

display(plot(η, v, xlabel = "η", ylabel = "v (m/s)", legend = false))
plot(x, v, xlabel = "Distance (m)", ylabel = "v (m/s)", legend = false)

# solveflat(v₀, vf, mytrack, myres, ρ, T)