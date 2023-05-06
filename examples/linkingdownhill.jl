using DifferentialEquations
using Plots
using OptimalTrainControl
using Roots
import ForwardDiff: derivative

function odefun!(du, u, p, x)
    # t, v, μ₂ = u
    t, v, η = u

    # η = μ₂ / v - 1.
    # ζ = η + 1 - ρ

    # π₁ = η > 0. ? η : 0.
    # π₂ = ζ ≥ 0. ? 0. : -ζ

    # @show v, μ₂, x, mycontrol(v, μ₂, x)

    du[1] = inv(v)
    du[2] = (p.u(u, p, x) - p.r(u, p, x) + 
        p.g(u, p, x)) * inv(v)
    # du[3] = -ψ(myresistance, V) / v^2 + μ₂ * du[2] / v + 
    #     μ₂ * derivative(x -> resistance(myresistance, x), v) / v - 
    #     π₁ * derivative(x -> 1/max(5,x), v) 
    #     π₂ * derivative(x -> -1/max(5,x), v)
    if η ≥ p.ρ - 1
        du[3] = (ψ(myresistance, v) - v^2 * derivative(v -> p.u([t, v, η], p, x), v)) * η / v^3 + (ψ(myresistance,v) - ψ(myresistance,V))/v^3
    else
        du[3] = (ψ(myresistance, v) - v^2 * derivative(v -> p.u([t, v, η], p, x), v)) * η / v^3 + (ψ(myresistance,v) - ψ(myresistance,V))/v^3 - (1 - p.ρ) * derivative(v -> p.u([t, v, η], p, x), v) / v
    end
end

condition_lowspeed(u, t, int) = u[2] - 1e-2
affect_lowspeed!(int) = terminate!(int)

function condition_modeswitch(u, t, int)
    if int.p.currentmode == :MaxP
        u[3]
    elseif int.p.currentmode == :Coast
        u[3] - (int.p.ρ - 1)
    elseif int.p.currentmode == :MaxB
        u[3] - (int.p.ρ - 1)
    end
end
function affect_modeswitch!(int)
    if int.p.currentmode == :Coast
        int.p.currentmode = :MaxP
    elseif int.p.currentmode == :MaxB
        int.p.currentmode = :Coast
    end
end
function affect_neg_modeswitch!(int)
    if int.p.currentmode == :MaxP
        int.p.currentmode = :Coast
    elseif int.p.currentmode == :Coast
        int.p.currentmode = :MaxB
    end
end

function mycontrol(u, p, x)
    if p.currentmode == :MaxP
        3. / u[2]
    elseif p.currentmode == :Coast
        0
    elseif p.currentmode == :MaxB
        - 3. / max(5., u[2])
    end
end

function condition(u, t, int)
    # @show t, u[3]
    u[3]
end

function solve_regular(u0, span, p0, seg2)
    cb_modeswitch = ContinuousCallback(condition_modeswitch, affect_modeswitch!, affect_neg_modeswitch!)
    cb_lowspeed = ContinuousCallback(condition_lowspeed, affect_lowspeed!)

    tstops = [r[:Distance] for r in eachrow(steephilltrack.waypoints)]
    
    cb = ContinuousCallback(condition, int -> int.t < seg2.start ? nothing : terminate!(int))
    cbS = CallbackSet(cb, cb_modeswitch, cb_lowspeed)
    
    prob = ODEProblem(odefun!, u0, span, p0)
    sol = solve(prob;alg_hints = [:stiff], tstops, callback = cbS, d_discontinuities = tstops)
end

function try_link(x0, seg2, initmode)
    p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, initmode)
    sol = solve_regular([0.0, V, 0.0], (x0, seg2.finish), p0, seg2)    

    v = sol[2,:]
    η = sol[3,:]
    x = sol.t

    # @show x[end]

    if x[end] == seg2.finish
        sign(η[end])*Inf
    elseif v[end] ≈ 1e-2
        -Inf
    else
        v[end] - V
    end
end

f(xprev, xnow, yprev, γ) = tan(asin(γ / -9.81))*(xnow - xprev) + yprev

function getys(γs, X)
    @assert length(X) == length(γs) + 1 "Lengths don't match!"
    ret = []
    for idx in eachindex(X)
        if idx == 1
            push!(ret, 0.0)
        else
            push!(ret, f(X[idx-1],X[idx],ret[idx-1],γs[idx-1]))
        end
    end
    return ret
end

trackX = [0.0,5e3,5.1e3, 10e3]
trackY = [0.0, 0.0, -10.0, -10.0]

steephilltrack = HillyTrack(trackX, Vector{Float64}(trackY))

myresistance = DavisResistance(6.75e-3, 0., 5e-5)


V = 20.0
segs = getmildsegments(steephilltrack, V, myresistance, x -> 3/x)

ρ = 0
# x0 = segs[2].finish - 383z
# x0 = segs[2].start
# xf = segs[2].finish

# vals = [try_link(x, segs[3]) for x in x0:5:xf]
# scatter(x0:5:xf,vals)
# μ₂0 = v₀ * ( (E(myresistance, V, v₀) - E(myresistance, V, V)) / 
        # (1. /max(5., v₀) - resistance(myresistance, v₀) + getgradientacceleration(steephilltrack, x0)) + 1.)

# display(plot(x, v, title = "Speed (m/s)"))

# plot(x, η, title = "η")
# plot!(twinx(), steephilltrack; alpha = 0.5)
p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, :Coast)
sol = solve_regular([0.0,V,0.0], (segs[2].finish-1e3, segs[3].finish), p0, segs[3])

plot(sol.t, sol[3,:])

# xs = segs[2].start:100:segs[2].finish

# vals = [clamp(try_link(x, segs[3], :Coast), -V, +V) for x in xs]
# scatter(xs, vals)

x_opt = find_zero(x -> try_link(x, segs[3], :Coast), [segs[2].start,segs[2].finish])

sol_opt = solve_regular([0.0,V,0.0],(x_opt,segs[3].finish),p0,segs[3])

plot(sol_opt.t,sol_opt[2,:])