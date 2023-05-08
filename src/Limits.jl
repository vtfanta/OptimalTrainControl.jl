@reexport using OptimalTrainControl

export linkunderspeedlimit

function linkunderspeedlimit(seg1::Segment, seg2::Segment, modelparams::ModelParams)
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
            u_max(u[2])
        elseif p.currentmode == :Coast
            0
        elseif p.currentmode == :MaxB
            u_min(u[2])
        end
    end
    
    function solve_regular!(u0, span, p0, seg2, x0 = nothing)
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
                push!(switching_points, (int.t, int.p.currentmode))
            elseif idx == 2 && int.p.currentmode == :MaxB
                int.p.currentmode = :Coast
                push!(switching_points, (int.t, int.p.currentmode))
            end
        end
        function affect_neg_modeswitch!(int, idx)
            if idx == 1 && int.t ≥ seg2.start
                terminate!(int)
            end
    
            if idx == 1 && int.p.currentmode == :MaxP
                int.p.currentmode = :Coast
                push!(switching_points, (int.t, int.p.currentmode))
            elseif idx == 2 && int.p.currentmode == :Coast
                int.p.currentmode = :MaxB
                push!(switching_points, (int.t, int.p.currentmode))
            end
        end
    
        function affect_speedlimit!(int)
            int.u[3] += x0
            if int.u[3] > 0
                int.p.currentmode = :MaxP
            elseif p0.ρ - 1 ≤ int.u[3] ≤ 0
                int.p.currentmode = :Coast
            else
                int.p.currentmode = :MaxB
            end
        end

        switching_points = []
    
        cb_modeswitch = VectorContinuousCallback(condition_modeswitch, affect_modeswitch!, affect_neg_modeswitch!, 2)
        cb_lowspeed = ContinuousCallback(condition_lowspeed, affect_lowspeed!)
    
        tstops = [r[:Distance] for r in eachrow(steephilltrack.waypoints)]
    
        if !isnothing(x0)
            cb_speedlimit = ContinuousCallback((u, t, int) -> u[2] - speedlimit(t), affect_speedlimit!)
            cbS = CallbackSet(cb_modeswitch, cb_lowspeed, cb_speedlimit)
        else
            cbS = CallbackSet(cb_modeswitch, cb_lowspeed)
        end
        
        prob = ODEProblem(odefun!, u0, span, p0)
        sol = solve(prob;alg_hints = [:stiff], tstops, callback = cbS, d_discontinuities = tstops, 
            dtmax = 10)
        return sol, switching_points
    end
    
    function try_link(x0, seg2, initmode, linkmode = :normal, start = 0)
        if linkmode == :normal
            p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
            (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, initmode)
            sol, _ = solve_regular!([0.0, V, 0.0], (x0, seg2.finish), p0, seg2)    
    
            v = sol[2,:]
            η = sol[3,:]
            x = sol.t
    
            @show x0
            # display(plot(x, v))
    
            if !checkspeedlimit(sol, speedlimit)
                return Inf
            end
    
            if x[end] ≈ seg2.finish
                sign(η[end]) * Inf
            elseif v[end] ≤ 0.1
                -Inf
            else
                if seg2.mode == :HoldP
                    v[end] - seg2.holdspeed
                elseif seg2.mode == :HoldR
                    v[end] - seg2.holdspeed
                end
            end
        elseif linkmode == :speedlimit
            p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
            (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, initmode)
            sol, _ = solve_regular!([0.0, V, 0.0], (start, seg2.finish), p0, seg2, x0)
    
            # @show sol.t[end]
            # @show sol[2,end],sol[3,end]
            # display(plot(sol.t, sol[3,:]))
            # display(vline!([seg2.start, seg2.finish]))
    
            if sol.t[end] ≈ seg2.finish
                sign(sol[3,end]) * Inf
            elseif sol[2,end] ≤ 0.1
                - Inf
            else
                sol[2,end] - seg2.holdspeed
            end
        end
    end

    steephilltrack = modelparams.track
    u_max = modelparams.umax
    u_min = modelparams.umin
    V = modelparams.V
    myresistance = modelparams.resistance
    vᵢ = modelparams.vᵢ
    vf = modelparams.vf
    ρ = modelparams.ρ
    speedlimit = modelparams.speedlimit
    @assert !isnothing(speedlimit)

    segs = getmildsegments(modelparams)

    if seg1.mode == :HoldP

        domain = (seg1.start, seg1.finish - 1)

        nudge = getgradientacceleration(steephilltrack, seg1.finish + 1) < 0 ? :MaxP : :Coast
    else # seg1.mode == :HoldR
        domain = (seg1.start, seg1.finish + 1)

        currg = getgradientacceleration(steephilltrack, seg1.finish - 1)
        nextg = getgradientacceleration(steephilltrack, seg1.finish + 1)

        nudge = (currg - nextg) > 0 ? :Coast : :MaxB
    end

    startingmode = nudge
    targetseg = seg2

    xopt = find_zero(x -> try_link(x, targetseg, startingmode), domain)
    p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
        (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, startingmode)
    sol, _ = solve_regular!([0.0,V,0.0], (xopt, targetseg.finish), p0, targetseg)
    # if maximum(sol[2,:]) ≈ speedlimit
        # println("CONDITION!")
    Δη = find_zero(η -> try_link(η, targetseg, startingmode, :speedlimit, xopt), (0, 5))
    p0 = OldModelParams(mycontrol, (u, p, x) -> resistance(myresistance, u[2]), 
    (u, p, x) -> getgradientacceleration(steephilltrack, x), ρ, startingmode)
    sol, points = solve_regular!([0.0,V,0.0], (xopt, targetseg.finish), p0, targetseg, Δη)
    # end
    pushfirst!(points, (xopt, nudge))
    push!(points, (sol.t[end], seg2.mode))
    return sol, points
end