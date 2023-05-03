using OptimalTrainControl
using Roots
using Plots
using DifferentialEquations
using NumericalIntegration
using ForwardDiff
using Optim

function playflatoptimal(m, r, len, T)

    MTscenario = MinimalTimeScenario(m, FlatTrack(len), 9.81, [1.0], [1.0])
    calculatecontrol!(MTscenario)
    MTsol = play(MTscenario)

    MTtotaltime = integrate(MTsol.t, inv.(MTsol[1,:]))

    @info "Time-optimal strategy reached the final point in $MTtotaltime secs."

    if MTtotaltime ≥ T
        error("The task is infeasible. Minimal time is $MTtotaltime s which is greater than the requested $T s.")
    end

    function f(V)
        @info "Simulating with V = $V m/s."
        _, t, _, _ = playflatfinishing(m, r, len, V[1], 0.0)
        t[end] - T
    end

    V = find_zero(f, (len / T / 2, 2len / T); xatol = 1e-2)

    x, t, v, η = playflatfinishing(m, r, len, V, 0.0)
end

function playflatfinishing(m::Model, r::Resistance, len::Real, V, lengthholdP = 0.0)

    x, t, v, η = playflat(m, r, len, V, lengthholdP)

    function f(LholdP)
        x, _, _, _ = playflat(m, r, len, V, LholdP)
        x[end] - len
    end

    if abs(x[end] - len) / len ≥ 1e-3
        newlengthholdP = find_zero((f,x->ForwardDiff.derivative(f, x)),
            len - x[end])
    else
        newlengthholdP = lengthholdP
    end

    @info "Found length of HoldP = $newlengthholdP m for V = $V m/s."

    x, t, v, η = playflat(m, r, len, V, newlengthholdP)
end

function playflat(m::Model, r::Resistance, len, V, lengthholdP)

    myresistance = r

    mymodel = m

    # For numerical purposes
    if lengthholdP ≤ typeof(lengthholdP)(1.0)
        holdPlength = typeof(lengthholdP)(1.0)
    else
        holdPlength = lengthholdP
    end

    myscenario = OptimalScenario(mymodel, FlatTrack(len), 9.81, typeof(V)[0.0, 1.0], [len, 1.0], V)

    # MAXIMUM POWER

    prob = ODEProblem(function (du, u, p, t) 
        du[1] = inv(u[2])
        du[2] = (mymodel.maxcontrol(u[2]) - resistance(myresistance, u[2])) * inv(u[2])
    end, myscenario.initialvalues, (0.0, len))

    condition1(u, t, int) = u[2] * one(V) - V
    affect1!(int) = terminate!(int)    
    cb = ContinuousCallback(condition1, affect1!)

    sol = solve(prob, callback = cb)
    t1 = sol[1,:]
    v1 = sol[2,:]

    E(Z, v) = ψ(myresistance, Z) / v + resistance(myresistance, v)

    rs = [resistance(myresistance, vk) for vk in v1]
    x1 = cumul_integrate(v1, v1 ./ (mymodel.maxcontrol.(v1) .- rs))
    ₊x₀ = x1[end]

    η1 = (E.(V, v1) .- E(V, V)) ./ (mymodel.maxcontrol.(v1) .- rs)

    # plot(η1, v1, color = :green, label = "MaxP")

    # CRUISING

    prob = ODEProblem(function (du, u, p, t) 
        du[1] = inv(u[2])
        du[2] = 0.0
    end, typeof(V)[t1[end], V], (x1[end], len))

    conditionC(u, t, int) = t*one(holdPlength) - (holdPlength + one(holdPlength)*x1[end])
    affectC!(int) = terminate!(int)
    cb = ContinuousCallback(conditionC, affectC!)

    sol = solve(prob, callback = cb)
    tC = sol[1,:]
    vC = sol[2,:]

    xC = x1[end] .+ V .* (tC .- tC[1])

    # COASTING
    prob = ODEProblem(function (du, u, p, t) 
        du[1] = inv(u[2])
        du[2] = -resistance(myresistance, u[2]) * inv(u[2])
    end, typeof(V)[tC[end], vC[end]], (xC[end], len))

    φ(v) = v * resistance(myresistance, v)
    Lᵥ(v) = ForwardDiff.derivative(φ, V) * (v - V) + φ(V)
    Uᵨᵥ = find_zero(x -> Lᵥ(x) - mymodel.ρ * φ(x), one(V) * 14.0)

    condition2(u, t, int) = u[2] * one(V) - Uᵨᵥ
    affect2!(int) = terminate!(int)    
    cb = ContinuousCallback(condition2, affect2!)

    sol = solve(prob, callback = cb)
    t2 = sol[1,:]
    v2 = sol[2,:]

    rs = [resistance(myresistance, vk) for vk in v2]

    x2 = xC[end] .+ cumul_integrate(v2, -v2 ./ rs)
    x₋ = x2[end]

    η2 = (E.(V, v2) .- E(V, V)) ./ (- rs)

    # plot!(η2, v2, color = :blue, label = "Coast")

    # MAXIMUM BRAKE

    prob = ODEProblem(function (du, u, p, t) 
        du[1] = inv(u[2])
        du[2] = (mymodel.mincontrol(u[2]) - resistance(myresistance, u[2])) * inv(u[2])
    end, typeof(V)[t2[end], v2[end]], (x2[end], len))

    condition3(u, t, int) = u[2] - 1e-4
    affect3!(int) = terminate!(int)    
    cb = ContinuousCallback(condition3, affect3!)

    sol = solve(prob, callback = cb)
    t3 = sol[1,:]
    v3 = sol[2,:]

    rs = [resistance(myresistance, vk) for vk in v3]

    x3 = x2[end] .+ cumul_integrate(v3, v3 ./ (mymodel.mincontrol.(v3) .- rs))

    η3 = mymodel.ρ * (E.(myscenario.W, v3) .- E.(myscenario.W, Uᵨᵥ)) ./ (mymodel.mincontrol.(v3) .- rs) .+ mymodel.ρ .- 1

    # plot!(η3, v3, color = :red, label = "MaxB")

    # plot([x1; x2; x3], [v1; v2; v3])

    return [x1; xC; x2; x3], [t1; tC; t2; t3], [v1; vC; v2; v3], [η1; zeros(length(vC)); η2; η3]
end

myresistance = DavisResistance(1e-2, 0., 1.5e-5)

mymodel = AlbrechtModel(myresistance, v -> inv(max(5., v)), v -> -inv(max(5., v)), 1e3, 0.5)

T = 1800

X = 30e3

V = 25.0

x, t, v, η = playflatoptimal(mymodel, myresistance, X, T)

plot(x, v)

# begin
#     Vs = [5k for k = 3.0:0.3:7.0]
#     ts = [playflatfinishing(mymodel, myresistance, X, Vs[k])[2][end] for k = eachindex(Vs)]
#     plot(Vs, ts .- T)
# end

# begin
# V = 26.38 # CORRECT FOR 30E3 M AND 1500 S
# # V = 26.
# x, t, v, η = playflatfinishing(mymodel, myresistance, X, V, 0.0)
# @show t[end]
# plot(x, v)
# end
