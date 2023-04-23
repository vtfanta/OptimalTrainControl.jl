@reexport using RailDynamics

struct NewModelParams
    umax
    umin
    resistance
    ρ
    track
    V
    vᵢ
    vf
end

mutable struct SolverParams
    modelparams::NewModelParams
    currentmode
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
        du[3] = (ψ(p.modelparams.resistance, v) - v^2 * derivative(v -> control(v, x, p), v)) * η / v^3 + (ψ(p.modelparams.resistance,v) - ψ(p.modelparams.resistance,V))/v^3 - (1 - ρ) * derivative(v -> control(v, x, p), v, v) / v
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

    if modelparams.ρ > 0
        res = modelparams.resistance
        V = modelparams.V
        W = find_zero(W -> modelparams.ρ * ψ(res, W) - ψ(res, V), (V, Inf))
    end

    # Link first segment to final segment
    if across
        sol, _ = solve_regular!([0.0, vinit, 0.0], (x0, finish(track)), p0, seg2)
        # display(plot(sol.t,sol[2,:]))
        if sol.retcode == ReturnCode.Terminated
            return -Inf
        else
            return sol[2,end] - modelparams.vf
        end
    end

    # Link first segment
    if isinf(seg2.start)
        sol, _ = solve_regular!([0.0, V, 0.0], (x0, seg2.finish), p0, seg2, false)    
        if sol.retcode == ReturnCode.Terminated # came to stop before finish
            return -Inf
        elseif sol.retcode == ReturnCode.Success
            return sol[2,end] - modelparams.vᵢ
        end
    end

    # Link final segment
    if isinf(seg2.finish)
        sol, _ = solve_regular!([0.0, V, 0.0], (x0, seg2.start), p0, seg2)
        if sol.retcode == ReturnCode.Terminated # came to stop before finish
            return -Inf
        elseif sol.retcode == ReturnCode.Success
            return sol[2,end] - modelparams.vf
        end        
    end

    sol, _ = solve_regular!([0.0, V, 0.0], (x0, seg2.finish), p0, seg2)    

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

        valleft = try_link(domain[1], seg2, nudge)
        valright = try_link(domain[2], seg2, nudge)

        if sign(valleft) == sign(valright)
            return nothing
        end

        xopt = find_zero(x -> try_link(x, seg2, nudge), domain; xatol = 1e-3)

        p = SolverParams(modelparams, nudge)

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
        p = SolverParams(modelparams, nudge)
        brakesol, pts = solve_regular!(u0, span, p, seg2)
        
        totalsolu = hcat(powersol[:,powersol.t .< switchingpoint], brakesol[:,:])
        totalsolt = vcat(powersol.t[powersol.t .< switchingpoint], brakesol.t)

        retsol = DiffEqBase.build_solution(ODEProblem(odefun!, [0.0, vᵢ, 1.0], (start(track), finish(track)), powerp), AutoTsit5(Rosenbrock23()), totalsolt, totalsolu, retcode = ReturnCode.Success)
        
        retpoints = vcat((start(track), :MaxP), (switchingpoint, :Coast), pts)

        return nothing
        return retsol, retpoints
    end
end