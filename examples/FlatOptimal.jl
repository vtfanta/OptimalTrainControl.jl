using RailDynamics
using Roots
using Plots
using DifferentialEquations
using NumericalIntegration
using ForwardDiff

myresistance = DavisResistance(1e-2, 0., 1.5e-5)

mymodel = AlbrechtModel(myresistance, v -> inv(max(5., v)), v -> -inv(max(5., v)), 1e3, 0.5)

V = 25.0

myscenario = OptimalScenario(mymodel, FlatTrack(30e3), 9.81, [0.0, 1.0], [30e3, 0.5], V)


# MAXIMUM POWER

prob = ODEProblem(function (du, u, p, t) 
    du[1] = inv(u[2])
    du[2] = (mymodel.maxcontrol(u[2]) - resistance(myresistance, u[2])) * inv(u[2])
end, myscenario.initialvalues, (0.0, 30e3+10))

condition(u, t, int) = u[2] - V
affect!(int) = terminate!(int)    
cb = ContinuousCallback(condition, affect!)

sol = solve(prob, callback = cb)
t1 = sol[1,:]
v1 = sol[2,:]


E(Z, v) = ψ(myresistance, Z) / v + resistance(myresistance, v)

rs = [resistance(myresistance, vk) for vk in v1]
x1 = cumul_integrate(v1, v1 ./ (mymodel.maxcontrol.(v1) .- rs))

η1 = (E.(V, v1) .- E(V, V)) ./ (mymodel.maxcontrol.(v1) .- rs)

plot(η1, v1, color = :green, label = "MaxP")

# COASTING
prob = ODEProblem(function (du, u, p, t) 
    du[1] = inv(u[2])
    du[2] = -resistance(myresistance, u[2]) * inv(u[2])
end, [t1[end], v1[end]], (x1[end], 30e3+10))

φ(v) = v * resistance(myresistance, v)
Lᵥ(v) = ForwardDiff.derivative(φ, V) * (v - V) + φ(V)
Uᵨᵥ = find_zero(x -> Lᵥ(x) - mymodel.ρ * φ(x), 14.0)

condition(u, t, int) = u[2] - Uᵨᵥ
affect!(int) = terminate!(int)    
cb = ContinuousCallback(condition, affect!)

sol = solve(prob, callback = cb)
t2 = sol[1,:]
v2 = sol[2,:]

rs = [resistance(myresistance, vk) for vk in v2]

x2 = x1[end] .+ cumul_integrate(v2, -v2 ./ rs)

η2 = (E.(V, v2) .- E(V, V)) ./ (- rs)

plot!(η2, v2, color = :blue, label = "Coast")

# MAXIMUM BRAKE

prob = ODEProblem(function (du, u, p, t) 
    du[1] = inv(u[2])
    du[2] = (mymodel.mincontrol(u[2]) - resistance(myresistance, u[2])) * inv(u[2])
end, [t2[end], v2[end]], (x2[end], 30e3+10))

condition(u, t, int) = u[2] - 1.0
affect!(int) = terminate!(int)    
cb = ContinuousCallback(condition, affect!)

sol = solve(prob, callback = cb)
t3 = sol[1,:]
v3 = sol[2,:]

rs = [resistance(myresistance, vk) for vk in v3]

x3 = x2[end] .+ cumul_integrate(v3, v3 ./ (mymodel.mincontrol.(v3) .- rs))

η3 = mymodel.ρ * (E.(myscenario.W, v3) .- E.(myscenario.W, Uᵨᵥ)) ./ (mymodel.mincontrol.(v3) .- rs) .+ mymodel.ρ .- 1

plot!(η3, v3, color = :red, label = "MaxB")

# plot([x1; x2; x3], [v1; v2; v3])