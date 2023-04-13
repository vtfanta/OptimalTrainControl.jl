using DifferentialEquations
using DiffEqCallbacks
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
        # 1. / max(u[2],5)
        0.125
    elseif p.currentmode == :Coast
        0
    elseif p.currentmode == :MaxB
        # - 3. / u[2]
        -0.25
    end
end

function solve_regular!(u0, span, p0, seg2)
    condition_lowspeed(u, t, int) = u[2] - 1e-2
    affect_lowspeed!(int) = terminate!(int)
    function condition_modeswitch(out, u, t, int)
        out[1] = u[3]
        out[2] = u[3] - (int.p.ρ - 1)
        out[3] = u[2] - (seg2.mode == :HoldP ? V : "TODO")
    end
    function affect_modeswitch!(int, idx)
        if idx == 3 && int.t ≥ seg2.start
            terminate!(int)
            return
        end

        if idx == 1 && int.p.currentmode == :Coast
            int.p.currentmode = :MaxP
            push!(switching_points, (int.t, int.p.currentmode))
        elseif idx == 2 && int.p.currentmode == :MaxB
            int.p.currentmode = :Coast
            push!(switching_points, (int.t, int.p.currentmode))
        end
    end
    function affect_neg_modeswitch!(int, idx)
        if idx == 3 && int.t ≥ seg2.start
            terminate!(int)
            return
        end

        if idx == 1 && int.p.currentmode == :MaxP
            int.p.currentmode = :Coast
            push!(switching_points, (int.t, int.p.currentmode))
        elseif idx == 2 && int.p.currentmode == :Coast
            int.p.currentmode = :MaxB
            push!(switching_points, (int.t, int.p.currentmode))
        end
    end

    switching_points = []

    cb_modeswitch = VectorContinuousCallback(condition_modeswitch, affect_modeswitch!, affect_neg_modeswitch!, 3)
    cb_lowspeed = ContinuousCallback(condition_lowspeed, affect_lowspeed!)

    tstops = [r[:Distance] for r in eachrow(steephilltrack.waypoints)]

    cbS = CallbackSet(cb_modeswitch, cb_lowspeed)
    
    prob = ODEProblem(odefun!, u0, span, p0)
    sol = solve(prob;alg_hints = [:stiff], tstops, callback = cbS, d_discontinuities = tstops, 
        dtmax = 10)
    return sol, switching_points
end

function try_link(x0, seg2, initmode, across = false)
    p0 = ModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, initmode)

    # Link first segment to final segment
    if across
        # sol, _ = solve_regular!([0.0, vᵢ, x0], (start(steephilltrack), finish(steephilltrack)), p0, seg2)
        # if sol.retcode == ReturnCode.Terminated
        #     return -Inf
        # else
        #     return sol[2,end] - vf
        # end
        sol, _ = solve_regular!([0.0, vf, x0], (finish(steephilltrack), start(steephilltrack)), p0, seg2)
        # display(plot(sol.t,sol[2,:]))
        if sol.retcode == ReturnCode.Terminated
            return -Inf
        else
            return sol[2,end] - vᵢ
        end
    end

    # Link first segment
    if isinf(seg2.start)
        sol, _ = solve_regular!([0.0, V, 0.0], (x0, seg2.finish), p0, seg2)    
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
        sol, _ = solve_regular!([0.0, V, 0.0], (x0, seg2.start), p0, seg2)
        if sol.retcode == ReturnCode.Terminated
            # display(plot(sol.t, sol[2,:]))
            return -Inf
        elseif sol.retcode == ReturnCode.Success
            # display(plot(sol.t, sol[2,:]))
            return sol[2,end] - vf
        end        
    end

    sol, _ = solve_regular!([0.0, V, 0.0], (x0, seg2.finish), p0, seg2)    

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

        sol, points = solve_regular!(u0, (xopt, seg2.finish), p, seg2)
        pushfirst!(points, (xopt, nudge))
        push!(points, (sol.t[end],seg2.mode))
        @show points
        return sol, points

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

        sol, points = solve_regular!(u0, (xopt, seg2.start), p, seg1)
        pushfirst!(points, (seg2.start, nudge))
        push!(points, (xopt, seg2.mode))
        @show points
        return sol, points

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

        sol, points = solve_regular!(u0, (xopt, seg2.start), p, seg2)
        pushfirst!(points, (xopt, nudge))
        push!(points, (sol.t[end], nothing))
        @show points
        return sol, points

    elseif isinf(seg1.start) && isinf(seg2.finish) # link initial to final segment
        nudge = :MaxB
        η = find_zero(η0 -> try_link(η0, seg2, nudge, true), (-1-1e-5,-100))
        p = ModelParams(u, (u, _, _) -> resistance(res, u[2]),
        (_, _, x) -> getgradientacceleration(track, x), ρ, nudge)
        sol, points = solve_regular!([0.0, vf, η], (finish(track), start(track)), p, seg2)

        reverse!(points)
        correct_points = []
        for point in points
            if point[2] == :MaxP
                push!(correct_points, (point[1], :Coast))
            elseif point[2] == :Coast
                if sol(point[1]+1)[3] < 0
                    push!(correct_points, (point[1], :MaxB))
                else
                    push!(correct_points, (point[1], :MaxP))
                end
            elseif point[2] == :MaxB
                push!(correct_points, (point[1], :Coast))
            end
        end
        pushfirst!(correct_points, (start(track), :MaxP))
        push!(correct_points, (finish(track), nothing))
        @show correct_points
        return sol, correct_points
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

# trackX = [-10e3,-0.8e3,-0.2e3,0,1e3,1.8e3,1.9e3,3.1e3,4.2e3,4.6e3,10e3]
# trackY = getys([0,0.03,0,0,-0.03,0,0.03,0,-0.04,0], trackX)

# steephilltrack = HillyTrack(trackX, Vector{Float64}(trackY))

# myresistance = DavisResistance(1e-2,0,1.5e-5)

# V = 25.0
# vᵢ = 5.0
# vf = 3.0
# segs = getmildsegments(steephilltrack, V, myresistance, x -> 1/max(x,5))
# @show segs
# ρ = 0

trackX = [0,16e3,20e3,24e3,25e3,28e3,31e3,40e3]
trackY = [0,0,400,160,160,460,280,280]/9.81

steephilltrack = HillyTrack(trackX, trackY)

myresistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)

V = sqrt(2 * 63.27)
vᵢ = sqrt(2 * 2)
vf = vᵢ
ρ = 0

segs = getmildsegments(steephilltrack, V, myresistance, x -> 0.125)

# sol0, _ = link(segs[3], segs[4], steephilltrack, mycontrol, myresistance, ρ, V)
# display(plot(sol0.t, sol0[2,:]; color = phasecolor.(sol0[3,:], ρ), lw = 3, label = false,
#     ylabel = "Speed (m/s)"))
# display(plot!(twinx(), steephilltrack; alpha = 0.5, size = (1600/2,900/2)))

# sol1, _ = link(segs[2], segs[4], steephilltrack, mycontrol, myresistance, ρ, V)
# display(plot!(sol1.t, sol1[2,:]; color = phasecolor.(sol1[3,:], ρ), lw = 3, label = false))

# sol2, _ = link(segs[4], segs[5], steephilltrack, mycontrol, myresistance, ρ, V)
# display(plot!(sol2.t, sol2[2,:]; color = phasecolor.(sol2[3,:], ρ), lw = 3, label = false))

# sol3, _ = link(segs[5], segs[6], steephilltrack, mycontrol, myresistance, ρ, V)
# display(plot!(sol3.t, sol3[2,:]; color = phasecolor.(sol3[3,:], ρ), lw = 3, label = false))

# sol4, _ = link(segs[2], segs[7], steephilltrack, mycontrol, myresistance, ρ, V)
# display(plot!(sol4.t, sol4[2,:]; color = phasecolor.(sol4[3,:], ρ), lw = 3, label = false))

# sol5, _ = link(segs[1], segs[7], steephilltrack, mycontrol, myresistance, ρ, V)
# display(plot!(sol5.t, sol5[2,:]; color = phasecolor.(sol5[3,:], ρ), lw = 3, label = false))

N = length(segs)
links = Matrix{Any}(nothing, N, N)
linkages = Dict()

links[1, 2] = link(segs[1], segs[2], steephilltrack, mycontrol, myresistance, ρ, V)
if !isnothing(links[1,2])
    linkages[segs[2]] = [segs[1]]
end

for k = 2:N-1
    z = k
    while !isnothing(z) && z > 0
        candidatelink = link(segs[z], segs[k+1], steephilltrack, mycontrol, myresistance, ρ, V)
        if isnothing(candidatelink)
            # no link possible
            z -= 1
        else # link found
            links[z, k+1] = candidatelink
            if isnothing(get(linkages, segs[k+1], nothing))
                linkages[segs[k+1]] = [segs[z]]
            elseif !in(segs[z], linkages[segs[k+1]])
                push!(linkages[segs[k+1]], segs[z])
            end
            if z == 1
                break
            end
            for i in z-1:1
                if isnothing(links[i,z])
                    continue
                end
                if links[i,z][2][end][1] < links[z,k+1][2][1][1] && 
                    !in(segs[i], linkages[segs[k+1]])
                    push!(linkages[segs[k+1]], segs[i])
                end
            end
            z = findlast([!in(seg, linkages[segs[k+1]]) for seg in segs[1:z-1]])
        end
    end
end





