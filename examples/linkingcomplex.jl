using DifferentialEquations
import NumericalIntegration: integrate
using Plots
using RailDynamics
using Roots
import ForwardDiff: derivative

function odefun!(du, u, p, x)
    t, v, η = u

    du[1] = inv(v)
    du[2] = (p.u(u, p, x) - p.r(u, p, x) + 
        p.g(u, p, x)) * inv(v)

    if p.currentmode != :MaxB
        du[3] = (ψ(myresistance, v) - v^2 * derivative(v -> p.u([t, v, η], p, x), v)) * η / v^3 + (ψ(myresistance,v) - ψ(myresistance,V))/v^3
    else
        du[3] = (ψ(myresistance, v) - v^2 * derivative(v -> p.u([t, v, η], p, x), v)) * η / v^3 + (ψ(myresistance,v) - ψ(myresistance,V))/v^3 - (1 - p.ρ) * derivative(v -> p.u([t, v, η], p, x), v) / v
    end
end

function mycontrol(u, p, x)
    if p.currentmode == :MaxP
        1. / max(u[2],5)
    elseif p.currentmode == :Coast
        0
    elseif p.currentmode == :MaxB
        - 3. / u[2]
    end
end

function solve_regular!(u0, span, p0, seg2)
    condition_lowspeed(u, t, int) = u[2] - 1e-2
    affect_lowspeed!(int) = terminate!(int)
    function condition_modeswitch(out, u, t, int)
        out[1] = u[3]
        out[2] = u[3] - (int.p.ρ - 1)
    end
    function affect_modeswitch!(int, idx)
        if idx == 1 && int.t ≥ seg2.start
            terminate!(int)
        end

        if idx == 1 && int.p.currentmode == :Coast
            int.p.currentmode = :MaxP
        elseif idx == 2 && int.p.currentmode == :MaxB
            int.p.currentmode = :Coast
        end
        # @show int.t, int.p.currentmode
    end
    function affect_neg_modeswitch!(int, idx)
        if idx == 1 && int.t ≥ seg2.start
            terminate!(int)
        end

        if idx == 1 && int.p.currentmode == :MaxP
            int.p.currentmode = :Coast
        elseif idx == 2 && int.p.currentmode == :Coast
            int.p.currentmode = :MaxB
        end
        # @show int.t, int.p.currentmode
    end

    cb_modeswitch = VectorContinuousCallback(condition_modeswitch, affect_modeswitch!, affect_neg_modeswitch!, 2)
    cb_lowspeed = ContinuousCallback(condition_lowspeed, affect_lowspeed!)

    tstops = [r[:Distance] for r in eachrow(steephilltrack.waypoints)]
    
    cbS = CallbackSet(cb_modeswitch, cb_lowspeed)
    
    prob = ODEProblem(odefun!, u0, span, p0)
    solve(prob;alg_hints = [:stiff], tstops, callback = cbS, d_discontinuities = tstops, 
        dtmax = 10)
end

function try_link(x0, seg2, initmode, accross = false)
    p0 = ModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, initmode)

    # Link first segment to final segment
    # if accross

    # end

    # Link first segment
    if isinf(seg2.start)
        sol = solve_regular!([0.0, V, 0.0], (x0, seg2.finish), p0, seg2)    
        if sol.retcode == ReturnCode.Terminated
            # display(plot(sol.t, sol[2,:]))
            return -Inf
        elseif sol.retcode == ReturnCode.Success
            # display(plot(sol.t, sol[2,:]))
            return sol[2,end] - vᵢ
        end
    end

    # Link final segment
    if isinf(seg2.finish)
        sol = solve_regular!([0.0, V, 0.0], (x0, seg2.start), p0, seg2)
        if sol.retcode == ReturnCode.Terminated
            # display(plot(sol.t, sol[2,:]))
            return -Inf
        elseif sol.retcode == ReturnCode.Success
            # display(plot(sol.t, sol[2,:]))
            return sol[2,end] - vf
        end        
    end

    sol = solve_regular!([0.0, V, 0.0], (x0, seg2.finish), p0, seg2)    

    v = sol[2,:]
    η = sol[3,:]
    x = sol.t

    # @show x[end]

    if x[end] ≈ seg2.finish
        sign(η[end]) * Inf
    elseif v[end] ≤ 0.1
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

function local_energy(u0, x0, p0, seg2)
    sol = solve_regular!(u0, (x0, seg2.finish), p0, seg2)
    v = sol[2,:]
    x = sol.t
    J = integrate(x, E.(myresistance, V, v) .- E(myresistance, V, V))
end

function link(seg1, seg2, track, u, res::DavisResistance, ρ, V)
    if ρ > 0
        W = find_zero(W -> ρ * ψ(res, W) - ψ(res, V), (V, Inf))
    end

    @info "Linking $seg1 \t+\t $seg2"

    if !isinf(seg1.start) && !isinf(seg2.finish) # link inner segments
        if seg1.finish > seg2.start
            error("Segment 1 has to precede segment 2.")
        end

        domain = (seg1.start, seg1.finish - 1)

        nudge = getgradientacceleration(track, seg1.finish + 1) < 0 ? :MaxP : :Coast

        valleft = try_link(domain[1], seg2, nudge)
        valright = try_link(domain[2], seg2, nudge)

        if sign(valleft) == sign(valright)
            return nothing
        end

        xopt = find_zero(x -> try_link(x, seg2, nudge), domain; xatol = 1e-3)

        p = ModelParams(u, (u, _, _) -> resistance(res, u[2]),
        (_, _, x) -> getgradientacceleration(track, x), ρ, nudge)

        if seg1.mode == :HoldP
            u0 = [0.0, V, 0.0]
        elseif seg1.mode == :HoldR
            u0 = [0.0, W, ρ - 1.0]
        end

        return solve_regular!(u0, (xopt, seg2.finish), p, seg2)

    elseif isinf(seg1.start) && !isinf(seg2.finish) # link initial segment to inner segment
        if seg1.finish > seg2.start
            error("Initial segment has no preceding segment.")
        end
        
        if seg2.mode == :HoldP
            nudge = vᵢ < V ? :MaxP : :Coast
        elseif seg2.mode == :HoldR
            nudge = vᵢ < W ? :Coast : :MaxB
        end

        domain = (seg2.start, seg2.finish)

        valleft = try_link(domain[1], seg1, nudge)
        valright = try_link(domain[2], seg1, nudge)

        if sign(valleft) == sign(valright)
            return nothing
        end

        xopt = find_zero(x -> try_link(x, seg1, nudge), domain; xatol = 1e-3)

        p = ModelParams(u, (u, _, _) -> resistance(res, u[2]),
        (_, _, x) -> getgradientacceleration(track, x), ρ, nudge)

        if seg2.mode == :HoldP
            u0 = [0.0, V, 0.0]
        elseif seg2.mode == :HoldR
            u0 = [0.0, W, ρ - 1.0]
        end

        return solve_regular!(u0, (xopt, seg2.start), p, seg1)

    elseif isinf(seg2.finish) && !isinf(seg1.start) # link final segment to inner segment
        if seg2.start < seg1.finish
            error("Final segment has no succeeding segment.")
        end

        if seg1.mode == :HoldP
            nudge = vf < V ? :Coast : :MaxP
        elseif seg1.mode == :HoldR
            nudge = vf < W ? :MaxB : :Coast
        end

        domain = (seg1.start, seg1.finish)

        valleft = try_link(domain[1], seg2, nudge)
        valright = try_link(domain[2], seg2, nudge)

        if sign(valleft) == sign(valright)
            return nothing
        end

        xopt = find_zero(x -> try_link(x, seg2, nudge), domain; xatol = 1e-3)

        p = ModelParams(u, (u, _, _) -> resistance(res, u[2]),
        (_, _, x) -> getgradientacceleration(track, x), ρ, nudge)

        if seg1.mode == :HoldP
            u0 = [0.0, V, 0.0]
        elseif seg1.mode == :HoldR
            u0 = [0.0, W, ρ - 1.0]
        end

        return solve_regular!(u0, (xopt, seg2.start), p, seg2)

    elseif isinf(seg1.start) && isinf(seg2.finish) # link initial to final segment

    end
end

function phasecolor(η, ρ)
    if η > 0
        return :green
    elseif ρ - 1 < η < 0
        return :grey
    else
        return :red
    end
end

trackX = [-10e3,-0.8e3,-0.2e3,0,1e3,1.8e3,1.9e3,3.1e3,4.2e3,4.6e3,10e3]
trackY = getys([0,0.03,0,0,-0.03,0,0.03,0,-0.04,0], trackX)

steephilltrack = HillyTrack(trackX, Vector{Float64}(trackY))

myresistance = DavisResistance(1e-2,0,1.5e-5)

V = 25.0
vᵢ = 3.0
vf = 2.0
segs = getmildsegments(steephilltrack, V, myresistance, x -> 1/max(x,5))
@show segs
ρ = 0

sol0 = link(segs[1], segs[2], steephilltrack, mycontrol, myresistance, ρ, V)
display(plot(sol0.t, sol0[2,:]; color = phasecolor.(sol0[3,:], ρ), lw = 3, label = false,
    ylabel = "Speed (m/s)"))
display(plot!(twinx(), steephilltrack; alpha = 0.5, size = (1600/2,900/2)))

sol1 = link(segs[2], segs[3], steephilltrack, mycontrol, myresistance, ρ, V)
display(plot!(sol1.t, sol1[2,:]; color = phasecolor.(sol1[3,:], ρ), lw = 3, label = false))

sol2 = link(segs[3], segs[5], steephilltrack, mycontrol, myresistance, ρ, V)
display(plot!(sol2.t, sol2[2,:]; color = phasecolor.(sol2[3,:], ρ), lw = 3, label = false))

sol3 = link(segs[5], segs[6], steephilltrack, mycontrol, myresistance, ρ, V)
display(plot!(sol3.t, sol3[2,:]; color = phasecolor.(sol3[3,:], ρ), lw = 3, label = false))

sol4 = link(segs[2], segs[7], steephilltrack, mycontrol, myresistance, ρ, V)
display(plot!(sol4.t, sol4[2,:]; color = phasecolor.(sol4[3,:], ρ), lw = 3, label = false))

# begin
# p0 = ModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
#     (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, :MaxP)
# sol = solve_regular!([0.0, vᵢ, 0.674740], (start(steephilltrack), finish(steephilltrack)), p0, segs[7])
# @show sol.t[end]
# plot(sol.t,sol[2,:])
# end