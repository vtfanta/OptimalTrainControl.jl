using RailDynamics
using DifferentialEquations
import ForwardDiff: derivative
import Roots: find_zero

flattrack = FlatTrack(20e3)

myresistance = DavisResistance(1e-2, 0., 1.5e-5)

# Cruising speed
V = 17.

# Proportion of recuperated braking energy
ρ = 0.5

# Initial speed
v₀ = 1.
# Final speed
vf = 1.

maxcontrol(v) = 1. / max(5., v)
mincontrol(v) = -1. / max(5., v)

function mycontrol(v, μ₂, x)
    if μ₂ > v
        maxcontrol(v)
    elseif μ₂ ≈ v
        resistance(myresistance, v) - getgradientacceleration(flattrack, x)
    elseif v > μ₂ > ρ * v
        0.
    elseif ρ * v > μ₂
        mincontrol(v)
    end
end

function odefun!(du, u, p, t)
    t = u[1]
    v = u[2]
    μ₂ = u[3]

    η = μ₂ / v - 1.
    ζ = η + 1 - ρ

    π₁ = η > 0. ? η : 0.
    π₂ = ζ ≥ 0. ? 0. : -ζ

    du[1] = inv(v)
    du[2] = (mycontrol(v, μ₂, t) - resistance(myresistance, v) + 
        getgradientacceleration(flattrack, t)) * inv(v)
    du[3] = -ψ(myresistance, V) / v^2 + μ₂ * du[2] / v + 
        μ₂ * derivative(x -> resistance(myresistance, x), v) / v - 
        π₁ * derivative(x -> mycontrol(x, Inf, start(flattrack)), v) + 
        π₂ * derivative(x -> mycontrol(x, -Inf, start(flattrack)), v)
end

# Calculate the initial value of μ₂
μ₂0 = v₀ * ( (E(myresistance, V, v₀) - E(myresistance, V, V)) / 
    (maxcontrol(v₀) - resistance(myresistance, v₀) + getgradientacceleration(flattrack, start(flattrack))) + 1.)

# Terminating at reaching the final speed
condition(u, t, int) = u[2] ≤ vf
affect!(int) = terminate!(int)
cb = DiscreteCallback(condition, affect!)
prob = ODEProblem(odefun!, [0., v₀, μ₂0], (start(flattrack), finish(flattrack)))

sol = solve(prob, alg_hints = [:stiff], callback = cb)
x = sol.t
t = sol[1,:]
v = sol[2,:]
μ₂ = sol[3,:]

η = μ₂ ./ v .- 1.

plot(η, v, xlabel = "η", ylabel = "v (m/s)", legend = false)