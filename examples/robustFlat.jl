using RailDynamics
using DifferentialEquations
import ForwardDiff: derivative
import Roots: find_zero
using Plots

flattrack = FlatTrack(20e3)

myresistance = DavisResistance(1e-2, 0., 1.5e-5)

# Cruising speed
V = 15.

# Proportion of recuperated braking energy
ρ = 0.5

# Initial speed
v₀ = 1.
# Final speed
vf = 1.

function __solveflat(v₀, vf, V, track::Track, res::Resistance, ρ, justswitchingpoints = false,
    maxcontrol = x -> 1. / max(5., x), mincontrol = x -> -1. / max(5., x), holdPend = 0.0)
    function mycontrol(v, μ₂, x)
        if μ₂ > v
            maxcontrol(v)
        elseif holdPend > 0.0 && x ≤ holdPend && v ≈ V
           resistance(res, v) - getgradientacceleration(flattrack, x)
        elseif v > μ₂ > ρ * v
            0.
        elseif ρ * v > μ₂
            mincontrol(v)
        else
            0.
        end
    end
    
    function odefun!(du, u, p, x)
        t = u[1]
        v = u[2]
        μ₂ = u[3]
    
        η = μ₂ / v - 1.
        ζ = η + 1 - ρ
    
        π₁ = η > 0. ? η : 0.
        π₂ = ζ ≥ 0. ? 0. : -ζ

        # @show v, μ₂, x, mycontrol(v, μ₂, x)
    
        du[1] = inv(v)
        du[2] = (mycontrol(v, μ₂, x) - resistance(res, v) + 
            getgradientacceleration(track, x)) * inv(v)
        du[3] = -ψ(res, V) / v^2 + μ₂ * du[2] / v + 
            μ₂ * derivative(x -> resistance(res, x), v) / v - 
            π₁ * derivative(x -> mycontrol(x, Inf, start(track)), v) + 
            π₂ * derivative(x -> mycontrol(x, -Inf, start(track)), v)
    end

    # Calculate the initial value of μ₂
    # Assuming that v₀ < V!!
    μ₂0 = v₀ * ( (E(res, V, v₀) - E(res, V, V)) / 
        (maxcontrol(v₀) - resistance(res, v₀) + getgradientacceleration(track, start(track))) + 1.)

    function condition(out, u, t, _)
        # Terminating at reaching the final speed
        out[1] = u[2] - vf
        # Reaching cruising speed
        out[2] = u[2] - u[3]
        # Reaching end of the HoldP phase
        out[3] = t - holdPend
    end
    function affect!(int, idx)
        if idx == 1
            terminate!(int)
        elseif idx == 2 && holdPend > 0.0
            @info "Hitting V, to singular HoldP"
            int.u[2] = V
            int.u[3] = V
        elseif idx == 2 && holdPend ≈ 0.0
            @info "Hitting V, to regular Coast"
            @show int.u[2], int.u[3]
            int.u[2] = V
            int.u[3] = V - 1e-3
        elseif idx == 3
            @info "Ending HoldP, to regular Coast"
            int.u[2] = V
            int.u[3] = V - 1e-3
        end
    end
    cb = VectorContinuousCallback(condition, affect!, 3)

    prob = ODEProblem(odefun!, [0., v₀, μ₂0], (start(track), finish(track)))

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

function _solveflat(v₀, vf, track, res, ρ, justswitchingpoints = false, 
    maxcontrol = x -> 1. / max(5., x), mincontrol = x -> -1. / max(5., x))
    
    x, t, v, η = __solveflat(v₀, vf, V, track, res, ρ, justswitchingpoints, maxcontrol, mincontrol)

    display(plot(x,v))

    holdPend = finish(track) - x[end]
    @show holdPend, x[end]
    
    x, t, v, η = __solveflat(v₀, vf, V, track, res, ρ, justswitchingpoints, maxcontrol, mincontrol, holdPend)
end

x, t, v, η = _solveflat(v₀, vf, flattrack, myresistance, ρ)

η_alg = (E.(myresistance, V, v) .- E(myresistance, V, V)) ./ (1.0 ./ max.(5., v) .- resistance.(myresistance, v))

plot(η, v, xlabel = "η", ylabel = "v (m/s)", legend = false)