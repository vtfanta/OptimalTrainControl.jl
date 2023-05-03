using DifferentialEquations
using DiffEqCallbacks
import NumericalIntegration: integrate
using Plots
using OptimalTrainControl
using Roots
import ForwardDiff: derivative

function odefun!(du, u, p, x)
    t, v, η = u

    du[1] = inv(v)
    du[2] = (p.u(u, p.currentmode, x) - p.r(u, p, x) + 
        p.g(u, p, x)) * inv(v)

    if p.currentmode != :MaxB
        du[3] = (ψ(myresistance, v) - v^2 * derivative(v -> p.u([t, v, η], p.currentmode, x), v)) * η / v^3 + (ψ(myresistance,v) - ψ(myresistance,V))/v^3
    else
        du[3] = (ψ(myresistance, v) - v^2 * derivative(v -> p.u([t, v, η], p.currentmode, x), v)) * η / v^3 + (ψ(myresistance,v) - ψ(myresistance,V))/v^3 - (1 - p.ρ) * derivative(v -> p.u([t, v, η], p.currentmode, x), v) / v
    end
end

function mycontrol(u, p, x)
    if p == :MaxP
        # 1. / max(u[2],5)
        0.125
    elseif p == :Coast
        0
    elseif p == :MaxB
        # - 3. / u[2]
        -0.25
    end
end

function solve_regular!(u0, span, p0, seg2, shouldswitch = true)
    condition_lowspeed(u, t, int) = u[2] - 1e-2
    affect_lowspeed!(int) = terminate!(int)
    function condition_modeswitch(out, u, t, int)
        out[1] = u[3]
        out[2] = u[3] - (int.p.ρ - 1)
        out[3] = u[2] - (seg2.mode == :HoldP ? V : "TODO")
    end
    function affect_modeswitch!(int, idx)
        if idx == 3 && int.t ≥ seg2.start && !isinf(seg2.start)
            terminate!(int)
            return
        end

        if shouldswitch
            if idx == 1 && int.p.currentmode == :Coast
                int.p.currentmode = :MaxP
                push!(switching_points, (int.t, int.p.currentmode))
            elseif idx == 2 && int.p.currentmode == :MaxB
                int.p.currentmode = :Coast
                push!(switching_points, (int.t, int.p.currentmode))
            end
        end
    end
    function affect_neg_modeswitch!(int, idx)
        if idx == 3 && int.t ≥ seg2.start && !isinf(seg2.start)
            terminate!(int)
            return
        end

        if shouldswitch
            if idx == 1 && int.p.currentmode == :MaxP
                int.p.currentmode = :Coast
                push!(switching_points, (int.t, int.p.currentmode))
            elseif idx == 2 && int.p.currentmode == :Coast
                int.p.currentmode = :MaxB
                push!(switching_points, (int.t, int.p.currentmode))
            end
        end
    end

    switching_points = []

    cb_lowspeed = ContinuousCallback(condition_lowspeed, affect_lowspeed!)

    tstops = [r[:Distance] for r in eachrow(steephilltrack.waypoints)]

    cb_modeswitch = VectorContinuousCallback(condition_modeswitch, affect_modeswitch!, affect_neg_modeswitch!, 3)
    cbS = CallbackSet(cb_modeswitch, cb_lowspeed)
    
    prob = ODEProblem(odefun!, u0, span, p0)
    sol = solve(prob;alg_hints = [:stiff], tstops, callback = cbS, d_discontinuities = tstops,)
        # dtmax = 100)
    return sol, switching_points
end

function try_link(x0, seg2, initmode, across = false, vinit = V)
    p0 = ModelParams((u,p,x) -> mycontrol(u,p,x), (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, initmode)

    # Link first segment to final segment
    if across
        sol, _ = solve_regular!([0.0, vinit, 0.0], (x0, finish(steephilltrack)), p0, seg2)
        # display(plot(sol.t,sol[2,:]))
        if sol.retcode == ReturnCode.Terminated
            return -Inf
        else
            return sol[2,end] - vf
        end
    end

    # Link first segment
    if isinf(seg2.start)
        sol, _ = solve_regular!([0.0, V, 0.0], (x0, seg2.finish), p0, seg2, false)    
        if sol.retcode == ReturnCode.Terminated
            # display(plot(sol.t, sol[2,:]))
            return -Inf
        elseif sol.retcode == ReturnCode.Success
            # display(plot(sol.t, sol[2,:]))
            @show sol[2,end] - vᵢ
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
        # display(plot(sol.t, sol[2,:]))
        sign(η[end]) * Inf
    elseif v[end] ≤ 0.1
        # display(plot(sol.t, sol[2,:]))
        -Inf
    else
        # display(plot(sol.t, sol[2,:]))
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
        # @show points

        if sol[2,end] ≤ 0.1
            # display(plot(sol.t, sol[2,:]))
            return nothing
        else
            return sol, points
        end

    elseif isinf(seg1.start) && !isinf(seg2.finish) # link initial segment to inner segment
        if seg1.finish > seg2.start
            error("Initial segment has no preceding segment.")
        end
        
        # if seg2.mode == :HoldP
        #     nudge = vᵢ < V ? :MaxP : :Coast
        # elseif seg2.mode == :HoldR
        #     nudge = vᵢ < W ? :Coast : :MaxB
        # end

        # domain = (seg2.start, seg2.finish)

        # valleft = try_link(domain[1], seg1, nudge)
        # valright = try_link(domain[2], seg1, nudge)

        # if sign(valleft) == sign(valright)
        #     return nothing
        # end

        # xopt = find_zero(x -> try_link(x, seg1, nudge), domain; xatol = 1e-3)

        # p = ModelParams(u, (u, _, _) -> resistance(res, u[2]),
        # (_, _, x) -> getgradientacceleration(track, x), ρ, nudge)

        # if seg2.mode == :HoldP
        #     u0 = [0.0, V, 0.0]
        # elseif seg2.mode == :HoldR
        #     u0 = [0.0, W, ρ - 1.0]
        # end

        # sol, points = solve_regular!(u0, (xopt, seg1.finish), p, seg1, false)
        # pushfirst!(points, (seg2.start, nudge))
        # push!(points, (xopt, seg2.mode))
        # @show points

        nudge = vᵢ < V ? :MaxP : :Coast
        u0 = [0.0, vᵢ, 1.0]
        p = ModelParams(u, (u, _, _) -> resistance(res, u[2]),
        (_, _, x) -> getgradientacceleration(track, x), ρ, nudge)
        sol, points = solve_regular!(u0, (start(track), seg2.finish), p, seg2, false)
        if sol[2,end] ≤ 0.1 || sol.retcode == ReturnCode.Success
            return nothing
        else
            return sol, [(start(track), nudge), (sol.t[end], seg2.mode)]
        end

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
        # @show points
        if sol[2,end] ≤ 0.1
            return nothing
        else
            return sol, points
        end

    elseif isinf(seg1.start) && isinf(seg2.finish) # link initial to final segment

        nudge = vᵢ < V ? :MaxP : :Coast
        powerp = ModelParams(u, (u, _, _) -> resistance(res, u[2]),
        (_, _, x) -> getgradientacceleration(track, x), ρ, nudge)
        powersol, _ = solve_regular!([0.0, vᵢ, 1.0], (start(track), finish(track)), 
            powerp, seg2, false)

        switchingpoint = find_zero(x -> try_link(x, seg2, :Coast, true, powersol(x)[2]), 
            (start(track), finish(track)))
       
        u0 = [powersol(switchingpoint)[1], powersol(switchingpoint)[2], 0.0]
        span = (switchingpoint, finish(track))
        nudge = :Coast
        p = ModelParams(u, (u, _, _) -> resistance(res, u[2]),
        (_, _, x) -> getgradientacceleration(track, x), ρ, nudge)
        brakesol, pts = solve_regular!(u0, span, p, seg2)
        
        totdistance = vcat(powersol.t[powersol.t .< switchingpoint], brakesol.t)
        tottime = powersol[1,powersol.t .< switchingpoint]
        totspeed = powersol[2,powersol.t .< switchingpoint]
        totη = powersol[3,powersol.t .< switchingpoint]
        append!(tottime, brakesol[1,:])
        append!(totspeed, brakesol[2,:])
        append!(totη, brakesol[3,:])
        totu = [[tottime[k], totspeed[k], totη[k]] for k in eachindex(totdistance)]

        retsol = DiffEqBase.build_solution(ODEProblem(odefun!, [0.0, vᵢ, 1.0], (start(track), finish(track)), powerp), AutoTsit5(Rosenbrock23()), totdistance, totu, retcode = ReturnCode.Success)
        
        retpoints = vcat((start(track), :MaxP), (switchingpoint, :Coast), pts)
        push!(retpoints, (finish(track), nothing))
        # display(plot(retsol.t, retsol.u[2,:]))

        return retsol, retpoints
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

function modecolor(xaxis, points)
    function mode2color(mode)
        if mode == :MaxP
            :green
        elseif mode == :Coast
            :grey
        elseif mode == :MaxB
            :red
        elseif mode == :HoldP
            :blue
        elseif mode == :HoldR
            :orange
        elseif isnothing(mode)
            :black
        end
    end
    colors = fill(:green, size(xaxis))
    pointidx = 1
    currcolor = mode2color(points[pointidx][2])
    for (xidx, x) in enumerate(xaxis)
        if pointidx < Base.length(points) && x ≥ points[pointidx+1][1]
            pointidx += 1
            currcolor = mode2color(points[pointidx][2])
        end
        colors[xidx] = currcolor
    end
    return colors
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

V = 6.5
# V = sqrt(2 * 63.27)
vᵢ = sqrt(2 * 2)
vf = vᵢ
ρ = 0

segs = getmildsegments(steephilltrack, V, myresistance, x -> 0.125, ρ)

function findchain(segs, track, control, res, ρ, V)
    function isinseg(x, seg)
        return seg.start ≤ x ≤ seg.finish
    end

    if ρ > 0
        W = find_zero(W -> ρ * ψ(res, W) - ψ(res, V), (V, Inf))
    end

    N = Base.length(segs)
    sols = Matrix{Any}(nothing, N, N)
    linkages = Dict()
    chains = Set()
    
    l = link(segs[1], segs[2], track, control, res, ρ, V)
    if !isnothing(l)
        linkages[2] = Set(1)
        union!(chains, [l[2]])
        sols[1,2] = l[1]
    end

    for k = 2:N-1
        z = k
        while z > 0
            # @show z, k+1
            l = link(segs[z], segs[k+1], track, control, res, ρ, V)
            if isnothing(l)
                # @show "NOPE"
                z -= 1
            else
                # @show "Link found"
                sols[z,k+1] = l[1]
                if z == 1 # Need to add new chain when linking to first segment
                    union!(chains, [l[2]])
                    linkages[k+1] = Set(z)
                end
                for c in chains
                    # @show chains
                    if segs[z].start ≤ c[end][1] ≤ segs[z].finish && c[end][1] < l[2][1][1]
                        newchain = vcat(c, l[2])
                        union!(chains, [newchain])
                        if (k+1) in keys(linkages)
                            union!(linkages[k+1], get(linkages, z, Set()))
                            union!(linkages[k+1], z)
                        else
                            linkages[k+1] = union(get(linkages, z, Set()), z)
                        end
                    else
                        linkages[k+1] = Set()
                    end
                end
                # Find furthest segment which is not linked to k+1
                z = findlast([!(s in linkages[k+1]) for s in eachindex(segs[1:z-1])])
                if isnothing(z)
                    break
                end
            end
        end
    end
    println(chains)
    for chain in chains
        if chain[1][1] == start(track) && chain[end][1] == finish(track) # Found overarching chain
            # Get indices of segments in the chain
            segsequence = filter(!isnothing,[c[2] == :HoldP || c[2] == :HoldR ? findfirst(isinseg.(c[1], segs)) : nothing for c in chain[2:end-1]])
            if isempty(segsequence)
                segsequence = [1, N]
            else
                pushfirst!(segsequence, 1)
                push!(segsequence, N)
            end

            @show segsequence

            # Build the complete solution
            K = Base.length(segsequence) # number of solution to combine
            if K > 2
                distances = [sols[segsequence[idx], segsequence[idx+1]].t for idx=1:K-1]
                times = [sols[segsequence[idx], segsequence[idx+1]][1,:] for idx=1:K-1]
                speeds = [sols[segsequence[idx], segsequence[idx+1]][2,:] for idx=1:K-1]
                ηs = [sols[segsequence[idx], segsequence[idx+1]][3,:] for idx=1:K-1]
                
                totdistance = vcat(distances...)
                tottime = times[1]
                totspeed = speeds[1]
                totη = ηs[1]
                for i = 2:K-1
                    timeoffset = (distances[i][1] - distances[i-1][end]) / 
                        (segs[i].mode == :HoldP ? V : W)
                    nexttime = times[i] .+ (tottime[end] + timeoffset)
                    nextspeed = speeds[i] 
                    nextη = ηs[i]

                    append!(tottime, nexttime)
                    append!(totspeed, nextspeed)
                    append!(totη, nextη)
                end

                # Construct the complete ODESolution
                totu = [[tottime[k], totspeed[k], totη[k]] for k in eachindex(totdistance)]
                prob = ODEProblem(odefun!, totu[1], (totdistance[1], totdistance[end]))
                sol = DiffEqBase.build_solution(prob, Tsit5(), totdistance, totu, 
                    retcode = ReturnCode.Success)

                return chain, sol
            else # just the initial and final segment linked
                return chain, sols[1,N]
            end
        end
    end
    print(chains)
    error("Chain from start to finish of the track not found.")
end

function calculatetime(chain)
    T_regular = sum([abs(c[1][1,end]) for c in chain])
    if length(chain) > 1
        T_singular = sum([chain[j][2][end][2] == :HoldP ? (chain[j+1][2][1][1] - chain[j][2][end][1]) / V : nothing for j = 1:length(chain)-1]) # TODO implement situation for :HoldR and W
    else
        T_singular = 0
    end
    T_regular + T_singular
end

function chain2sol(chain, track, control, res)
    @assert start(track) == chain[1][1]
    @assert finish(track) == chain[end][1]

    # parameter p serves as the current mode indicator
    # set p to the initial mode
    chainidx = 1
    p = chain[chainidx][2]

    switchingpoints = [c[1] for c in chain]

    function small_odefun!(du, u, p, x)
        t, v = u
        du[1] = inv(v)
        du[2] = (control(u, p, x) - resistance(res, u[2] + getgradientacceleration(track, x))) * 
            inv(v)
    end

    condition(u, x, int) = x ∈ switchingpoints
    function affect!(int)
        chainidx += 1
        if !isnothing(chain[chainidx][2])
            int.p = chain[chainidx][2]
        end
    end
    cb = DiscreteCallback(condition, affect!)

    trackchanges = [r[:Distance] for r in eachrow(track.waypoints)]
    tstops = sort(vcat(switchingpoints, trackchanges))

    prob = ODEProblem(small_odefun!, [0.0, vᵢ], (start(track), finish(track)), p)

    sol = solve(prob; alg_hints = [:stiff], tstops, d_discontinuities = tstops, callback = cb)

end

chain, sol = findchain(segs, steephilltrack, mycontrol, myresistance, ρ, V)
print(chain)
plot(sol.t, sol[2,:]; color = modecolor(sol.t, chain), lw = 2)