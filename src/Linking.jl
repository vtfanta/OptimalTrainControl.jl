@reexport using RailDynamics

export solve!, modecolor

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

function solve!(prob::TrainProblem; atol = 5)
    @unpack track, ρ, vᵢ, vf, T, resistance, umax, umin = prob
    if isa(track, FlatTrack)
        track = HillyTrack([0, track.X], [0,0])
    end

    function time_constraint(V)
        params = NewModelParams(umax, umin, resistance, ρ, track, V, vᵢ, vf)
        segs = getmildsegments(params)

        _, sol = findchain(segs, params)
        return sol[1,end] - T
    end

    # 140 m/s is above 500 km/h, so pretty much a guaranteed upper bound
    spanV = (RailDynamics.length(track) / T / 2, 140)

    Vopt = nothing
    try
        Vopt = find_zero(time_constraint, spanV; atol)
    catch e
        printstyled("The problem is likely infeasible. Try higher time of journey.";
            color = :red)
        return 
    end

    params = NewModelParams(umax, umin, resistance, ρ, track, Vopt, vᵢ, vf)
    chain, sol = findchain(getmildsegments(params), params)
    prob.states = sol
    prob.switchingpoints = chain
    return 
end

function control(v, x, params::SolverParams)
    mode = params.currentmode
    if mode == :MaxP
        params.modelparams.umax(v)
    elseif mode == :HoldP || mode == :HoldR
        resistance(params.modelparams.resistance, v) - getgradientacceleration(params.modelparams.track, x)
    elseif mode == :Coast
        0
    elseif mode == :MaxB
        params.modelparams.umin(v)
    end
end

function odefun!(du, u, p::SolverParams, x)
    t, v, η = u

    V = p.modelparams.V
    ρ = p.modelparams.ρ

    du[1] = inv(v)
    du[2] = (control(v, x, p) - resistance(p.modelparams.resistance, v) + 
        getgradientacceleration(p.modelparams.track, x)) * inv(v)

    if p.currentmode != :MaxB
        du[3] = (ψ(p.modelparams.resistance, v) - v^2 * derivative(v -> control(v, x, p), v)) * η / v^3 + (ψ(p.modelparams.resistance,v) - ψ(p.modelparams.resistance,V))/v^3
    else
        du[3] = (ψ(p.modelparams.resistance, v) - v^2 * derivative(v -> control(v, x, p), v)) * η / v^3 + (ψ(p.modelparams.resistance,v) - ψ(p.modelparams.resistance,V))/v^3 - (1 - ρ) * derivative(v -> control(v, x, p), v) / v
    end
end

function solve_regular!(u0, span, p::SolverParams, seg2::Segment, shouldswitch = true)
    ρ = p.modelparams.ρ
    res = p.modelparams.resistance
    V = p.modelparams.V
    track = p.modelparams.track

    if ρ > 0
        W = find_zero(W -> ρ * ψ(res, W) - ψ(res, V), (V, Inf))
    end

    condition_lowspeed(u, t, int) = u[2] - 1e-2
    affect_lowspeed!(int) = terminate!(int)
    function condition_modeswitch(out, u, t, int)
        out[1] = u[3]
        out[2] = u[3] - (ρ - 1)
        out[3] = u[2] - (seg2.mode == :HoldP ? V : W)
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

    tstops = [r[:Distance] for r in eachrow(track.waypoints)]

    cb_modeswitch = VectorContinuousCallback(condition_modeswitch, affect_modeswitch!, affect_neg_modeswitch!, 3)
    cbS = CallbackSet(cb_modeswitch, cb_lowspeed)
    
    prob = ODEProblem(odefun!, u0, span, p)
    sol = solve(prob;alg_hints = [:stiff], tstops, callback = cbS, d_discontinuities = tstops,)
    return sol, switching_points
end

function try_link(x0, seg2, initmode, modelparams::NewModelParams,  across = false, vinit = modelparams.V)
    p0 = SolverParams(modelparams, initmode)

    segs = getmildsegments(modelparams)
    seg1idx = findfirst([seg.start ≤ x0 ≤ seg.finish for seg in segs])
    if isnothing(seg1idx)
        seg1 = segs[1]
    else
        seg1 = segs[seg1idx]
    end

    V = modelparams.V
    if modelparams.ρ > 0
        res = modelparams.resistance
        W = find_zero(W -> modelparams.ρ * ψ(res, W) - ψ(res, V), (V, Inf))
    end

    # Link first segment to final segment
    if across
        sol, _ = solve_regular!([0.0, vinit, 0.0], (x0, finish(modelparams.track)), p0, seg2)
        if sol.retcode == ReturnCode.Terminated
            return -Inf
        else
            return sol[2,end] - modelparams.vf
        end
    end

    # Link first segment
    if isinf(seg2.start)
        u0 = seg2.mode == :HoldP ? [0.0, V, 0.0] : [0.0, W, modelparams.ρ - 1]
        sol, _ = solve_regular!(u0, (x0, seg2.finish), p0, seg2, false)    
        if sol.retcode == ReturnCode.Terminated # came to stop before finish
            return -Inf
        elseif sol.retcode == ReturnCode.Success
            return sol[2,end] - modelparams.vᵢ
        end
    end

    # Link final segment
    if isinf(seg2.finish)
        u0 = seg1.mode == :HoldP ? [0.0, V, 0.0] : [0.0, W, modelparams.ρ - 1]
        sol, _ = solve_regular!(u0, (x0, seg2.start), p0, seg2)
        if sol.retcode == ReturnCode.Terminated # came to stop before finish
            return -Inf
        elseif sol.retcode == ReturnCode.Success
            return sol[2,end] - modelparams.vf
        end        
    end

    # Link inner segments
    u0 = seg1.mode == :HoldP ? [0.0, V, 0.0] : [0.0, W, modelparams.ρ - 1]
    sol, _ = solve_regular!(u0, (x0, seg2.finish), p0, seg2)    

    v = sol[2,:]
    η = sol[3,:]
    x = sol.t

    if x[end] ≈ seg2.finish
        sign(η[end]) * Inf
    elseif v[end] ≤ 0.1
        -Inf
    else
        if seg2.mode == :HoldP
            v[end] - modelparams.V
        elseif seg2.mode == :HoldR
            v[end] - W
        end
    end
end

function link(seg1::Segment, seg2::Segment, modelparams::NewModelParams)
    V = modelparams.V
    res = modelparams.resistance
    ρ = modelparams.ρ
    track = modelparams.track
    vᵢ = modelparams.vᵢ
    vf = modelparams.vf

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

        valleft = try_link(domain[1], seg2, nudge, modelparams)
        valright = try_link(domain[2], seg2, nudge, modelparams)

        if sign(valleft) == sign(valright)
            return nothing
        end

        xopt = find_zero(x -> try_link(x, seg2, nudge, modelparams), domain; xatol = 1e-3)

        p = SolverParams(modelparams, nudge)

        if seg1.mode == :HoldP
            u0 = [0.0, V, 0.0]
        elseif seg1.mode == :HoldR
            u0 = [0.0, W, ρ - 1.0]
        end

        sol, points = solve_regular!(u0, (xopt, seg2.finish), p, seg2)
        pushfirst!(points, (xopt, nudge))
        push!(points, (sol.t[end],seg2.mode))

        if sol[2,end] ≤ 0.1
            return nothing
        else
            return sol, points
        end

    elseif isinf(seg1.start) && !isinf(seg2.finish) # link initial segment to inner segment
        if seg1.finish > seg2.start
            error("Initial segment has no preceding segment.")
        end
        
        nudge = vᵢ < V ? :MaxP : :Coast
        u0 = [0.0, vᵢ, 1.0]
        p = SolverParams(modelparams, nudge)
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

        valleft = try_link(domain[1], seg2, nudge, modelparams)
        valright = try_link(domain[2], seg2, nudge, modelparams)

        if sign(valleft) == sign(valright)
            return nothing
        end

        xopt = find_zero(x -> try_link(x, seg2, nudge, modelparams), domain; xatol = 1e-3)

        p = SolverParams(modelparams, nudge)

        if seg1.mode == :HoldP
            u0 = [0.0, V, 0.0]
        elseif seg1.mode == :HoldR
            u0 = [0.0, W, ρ - 1.0]
        end

        sol, points = solve_regular!(u0, (xopt, seg2.start), p, seg2)
        pushfirst!(points, (xopt, nudge))
        push!(points, (sol.t[end], nothing))
        if sol[2,end] ≤ 0.1
            return nothing
        else
            return sol, points
        end

    elseif isinf(seg1.start) && isinf(seg2.finish) # link initial to final segment
        nudge = vᵢ < V ? :MaxP : :Coast
        powerp = SolverParams(modelparams, nudge)
        powersol, _ = solve_regular!([0.0, vᵢ, 1.0], (start(track), finish(track)), 
            powerp, seg2, false)

        switchingpoint = find_zero(x -> try_link(x, seg2, :Coast,modelparams, 
            true, powersol(x)[2]), (start(track), finish(track)))
       
        u0 = [powersol(switchingpoint)[1], powersol(switchingpoint)[2], 0.0]
        span = (switchingpoint, finish(track))
        nudge = :Coast
        p = SolverParams(modelparams, nudge)
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
        return retsol, retpoints
    end
end

function findchain(segs, modelparams::NewModelParams)
    @show modelparams.V
    function isinseg(x, seg)
        return seg.start ≤ x ≤ seg.finish
    end

    ρ = modelparams.ρ
    V = modelparams.V
    track = modelparams.track
    res = modelparams.resistance

    if ρ > 0
        W = find_zero(W -> ρ * ψ(res, W) - ψ(res, V), (V, Inf))
    end

    N = Base.length(segs)
    sols = Matrix{Any}(nothing, N, N)
    linkages = Dict()
    chains = Set()
    
    l = link(segs[1], segs[2], modelparams)
    if !isnothing(l)
        linkages[2] = Set(1)
        union!(chains, [l[2]])
        sols[1,2] = l[1]
    end

# Solution can be connected to a chain if the chain precedes the solution and the endpoint of the 
# chain lies before the start of the solution.

    # for k = 2:N-1
    #     z = k
    #     while z > 0
    #         l = link(segs[z], segs[k+1], modelparams)
    #         if isnothing(l)
    #             z -= 1
    #         else
    #             sols[z, k+1] = l[1]
    #             if z == 1
    #                 union!(chains, [l[2]])
    #             end
    #             for c in chains
    #                 if segs[z].start ≤ c[end][1] ≤ segs[z].finish && c[end][1] < l[2][1][1]
    #                     newchain = vcat(c, l[2])
    #                     union!(chains, [newchain])
    #                     if (k+1) in keys(linkages)
    #                         union!(linkages[k+1], )
    #                     else

    #                     end
    #                 end
    #             end
    #         end
    #     end
    # end

    for k = 2:N-1
        z = k
        while z > 0
            # @show z, k+1
            l = link(segs[z], segs[k+1], modelparams)
            if isnothing(l)
                println("$z -> $(k+1) NOPE!")
                z -= 1
            else
                println("$z -> $(k+1) link found!")
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
                
                z = findlast([!(s in get(linkages, k+1, Set())) for s in eachindex(segs[1:z-1])])
                if isnothing(z)
                    break
                end
            end
        end
    end
    candidatechains = filter(c -> c[1][1] == start(track) && c[end][1] == finish(track), chains)
    chain = argmax(Base.length, candidatechains)
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
            # totη = ηs[1]
            for i = 2:K-1
                timeoffset = (distances[i][1] - distances[i-1][end]) / 
                    (segs[i].mode == :HoldP ? V : W)
                nexttime = times[i] .+ (tottime[end] + timeoffset)
                nextspeed = speeds[i] 
                # nextη = ηs[i]

                append!(tottime, nexttime)
                append!(totspeed, nextspeed)
                # append!(totη, nextη)
            end

            # Construct the complete ODESolution
            # totu = [[tottime[k], totspeed[k], totη[k]] for k in eachindex(totdistance)]
            totu = [[tottime[k], totspeed[k]] for k in eachindex(totdistance)]
            prob = ODEProblem(odefun!, totu[1], (totdistance[1], totdistance[end]))
            sol = DiffEqBase.build_solution(prob, Tsit5(), totdistance, totu, 
                retcode = ReturnCode.Success)

            return chain, sol
        else # just the initial and final segment linked
            totu = [[sols[1,N][1,k], sols[1,N][2,k]] for k in eachindex(sols[1,N].t)]
            totdistance = sols[1,N].t
            prob = ODEProblem(odefun!, totu[1], (totdistance[1], totdistance[end]))
            sol = DiffEqBase.build_solution(prob, Tsit5(), totdistance, totu, 
                retcode = ReturnCode.Success)
            return chain, sol
        end
    end
    # print(chains)
    error("Chain from start to finish of the track not found.")
end