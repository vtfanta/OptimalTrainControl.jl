using DifferentialEquations
using Plots
using OptimalTrainControl
using Roots

function _solve_flat(vᵢ, mytrack, myres, V, V′ = V)

    function maxcontrol(s, p, x)
        1 / max(5, s[2])
    end

    function mincontrol(s, p, x)
        -1 / max(5, s[2])
    end

    function mycontrol(s, p, x)
        if p.currentphase == :MaxP
            maxcontrol(s, p, x)
        elseif p.currentphase == :HoldP || p.currentphase == :HoldR
            p.r(s, p, x) - p.g(s, p, x)
        elseif p.currentphase == :Coast
            0
        elseif p.currentphase == :MaxB
            mincontrol(s, p, x)
        end
    end

    function condition(s, x, int)
        phase = int.p.currentphase

        # safeguard
        if s[2] ≤ 1e-2
            return 0
        end

        if phase == :MaxP
            s[2] - V′
        elseif phase == :Coast
            η̄ = (E(myres, V, s[2]) - E(myres, V, V′)) / -int.p.r(s, p, x)
            η̄ - (int.p.ρ - 1)
        else
            1
        end

    end

    function affect!(int)
        if int.u[2] ≤ 1e-2
            terminate!(int)
        end

        phase = int.p.currentphase

        if phase == :MaxP
            int.p.currentphase = :Coast
        elseif phase == :Coast
            int.p.currentphase = :MaxB
        end
    end

    p = OldModelParams(mycontrol, (s, p, x) -> resistance(myres, s[2]), 
        (s, p, x) -> 0, 0.5, :MaxP)

    s₀ = [0.0, vᵢ]

    W = find_zero(W -> p.ρ * ψ(myres, W) - ψ(myres, V), V′)
    Uᵨᵥ = find_zero(v -> φ′(myres, V) * (v - V) + φ(myres, V) - p.ρ * φ(myres, v), V)

    prob = ODEProblem(albrecht_odefun!, s₀, (start(mytrack), Inf), p)

    x = []
    t = []
    v = []
    η = []

    int = init(prob, alg_hints = [:stiff], 
            callback = ContinuousCallback(condition, affect!))

    for (s, x̄) in tuples(int)
        
        push!(t, s[1])
        push!(v, s[2])
        push!(x, x̄)

        phase = int.p.currentphase
        if phase == :MaxP || phase == :HoldP || phase == :Coast
            push!(η, (E(myres, V, s[2]) - E(myres, V, V′)) / 
                (int.p.u(s, int.p, x̄) - int.p.r(s, int.p, x̄)))
        elseif phase == :MaxB
            ζ = int.p.ρ * (E(myres, W, s[2]) - E(myres, W, Uᵨᵥ)) / 
                (int.p.u(s, int.p, x̄) - int.p.r(s, int.p, x̄))
            push!(η, ζ + int.p.ρ - 1)
        end
    end

    if abs(finish(mytrack) - x[end]) ≤ 5
        @info "Solution hit the end!"
    elseif V ≈ V′ && x[end] < finish(mytrack)
        @info "Inserting HoldP"
        holdPlength = finish(mytrack) - x[end]
        # @assert maximum(v) ≈ V "max(v) = $(maximum(v))"
        idx = findfirst(isapprox.(V′, v))
        insert!(t, idx, t[idx])
        t[idx+1:end] .+= holdPlength / V′
        insert!(v, idx, V′)
        insert!(x, idx, x[idx])
        x[idx+1:end] .+= holdPlength
        insert!(η, idx, 0.)
    elseif x[end] > finish(mytrack)
        @info "Passed the target distance"
    elseif x[end] < finish(mytrack)
        @info "Undershot the target distance"
    end
    x, t, v, η
end

function solve_flat(vᵢ, track, res, V)
    # V′ = V
    x, t, v, η = _solve_flat(vᵢ, track, res, V)

    # while x[end] > finish(track) + 5
    #     V′ *= 0.99
    #     @info "Overshoot. Reducing to V′ = $V′ m/s"
    #     x, t, v, η = _solve_flat(vᵢ, track, res, V, V′)
    # end
    if x[end] > finish(track) + 5
        V′ = find_zero(
                V′ -> _solve_flat(vᵢ, track, res, V, V′)[1][end] - finish(track), [2.,V])

        _solve_flat(vᵢ, track, res, V, V′)
    else
        x, t, v, η
    end
end

V = 50
x, t, v, η = solve_flat(1.0, FlatTrack(3e3), DavisResistance(1e-2, 0., 1.5e-5), V)

plot(x, v,
    label = false,
    xlabel = "Distance (m)",
    ylabel = "Speed (m/s)")