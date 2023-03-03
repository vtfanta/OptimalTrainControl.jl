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

function __solveflat(v₀, vf, V, track::Track, res::Resistance, ρ, justswitchingpoints = true,
    maxcontrol = x -> 1. / max(5., x), mincontrol = x -> -1. / max(5., x), holdPend = 0.0)
    function mycontrol(v, μ₂, x)
        if μ₂ > v
            maxcontrol(v)
        # elseif x ≤ holdPend
            # resistance(res, v) - getgradientacceleration(flattrack, x)
        elseif v > μ₂ > ρ * v
            0.
        elseif ρ * v > μ₂
            mincontrol(v)
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
    
        du[1] = inv(v)
        du[2] = (mycontrol(v, μ₂, t) - resistance(res, v) + 
            getgradientacceleration(track, t)) * inv(v)
        du[3] = -ψ(res, V) / v^2 + μ₂ * du[2] / v + 
            μ₂ * derivative(x -> resistance(res, x), v) / v - 
            π₁ * derivative(x -> mycontrol(x, Inf, start(track)), v) + 
            π₂ * derivative(x -> mycontrol(x, -Inf, start(track)), v)
    end

    # Calculate the initial value of μ₂
    # Assuming that v₀ < V!!
    μ₂0 = v₀ * ( (E(res, V, v₀) - E(res, V, V)) / 
        (maxcontrol(v₀) - resistance(res, v₀) + getgradientacceleration(track, start(track))) + 1.)

    # Terminating at reaching the final speed
    condition(u, _, _) = u[2] ≤ vf
    affect!(int) = terminate!(int)
    cb = DiscreteCallback(condition, affect!)
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

    holdPend = finish(track) - x[end]
    @show holdPend, x[end]
    
    x, t, v, η = __solveflat(v₀, vf, V, track, res, ρ, justswitchingpoints, maxcontrol, mincontrol, holdPend)
end

x, t, v, η = _solveflat(v₀, vf, flattrack, myresistance, ρ)

plot(η, v, xlabel = "η", ylabel = "v (m/s)", legend = false)