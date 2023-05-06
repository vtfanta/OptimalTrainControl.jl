### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 30410070-c901-11ed-3a7b-23fe70952345
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ 6f33c454-b47f-4435-9623-20b4783fea11
begin
	Pkg.add("PlutoUI")
end

# ╔═╡ 3b663ae4-d7ad-4e06-a3a8-88b54b52a290
using DifferentialEquations

# ╔═╡ d0e4a291-ff3b-484d-959a-d185a8d71073
using OptimalTrainControl

# ╔═╡ 425edfc2-c7d8-47c0-959f-11cd7bd58e57
using Roots

# ╔═╡ 544c1c3e-cfb5-43b8-b56b-fdc7c9d30c93
using Plots

# ╔═╡ 601ae6be-7809-4c62-851f-f549dab359d1
using PlutoUI

# ╔═╡ e7d708cc-5e2a-4abb-9172-8d4c045407c6
function __solve_flat(vᵢ, mytrack, myres, V, V′ = V)

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
        # @info "Solution hit the end!"
    elseif V ≈ V′ && x[end] < finish(mytrack)
        # @info "Inserting HoldP"
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
        # @info "Passed the target distance"
    elseif x[end] < finish(mytrack)
        # @info "Undershot the target distance"
    end
    x, t, v, η
end

# ╔═╡ e774edc5-4b9e-41ae-b86c-06f7b484f343
function _solve_flat(vᵢ, track, res, V)
    # V′ = V
    x, t, v, η = __solve_flat(vᵢ, track, res, V)

    # while x[end] > finish(track) + 5
    #     V′ *= 0.99
    #     @info "Overshoot. Reducing to V′ = $V′ m/s"
    #     x, t, v, η = _solve_flat(vᵢ, track, res, V, V′)
    # end
    if x[end] > finish(track) + 5
        V′ = find_zero(
                V′ -> __solve_flat(vᵢ, track, res, V, V′)[1][end] - finish(track), [2, V])

        __solve_flat(vᵢ, track, res, V, V′)
    else
        x, t, v, η
    end
end

# ╔═╡ a92c7259-91e1-4db9-840a-7f86aec18c1d
function solve_flat(vᵢ, track, res, T)
	V = find_zero(
		V -> _solve_flat(vᵢ, track, res, V)[2][end] - T, finish(track) / T
	)
	_solve_flat(vᵢ, track, res, V)
end

# ╔═╡ 5d6314fa-8e7c-4e22-8afe-11603be925bc
@bind V Slider(5:1:200)

# ╔═╡ 32dd4a90-39e7-440a-8802-bbbdd3062680
md"V = $V m/s"

# ╔═╡ 73445f47-8029-492c-927f-2da1f9b51a75
@bind V′ Slider(3:1:30)

# ╔═╡ 6a37fcbd-2ef5-49fd-8792-05c76e0720ed
md"V′ = $V′ m/s"

# ╔═╡ 1f9e3829-f18c-433c-8b7e-928bb7f92bd3
begin
	track = FlatTrack(5e3)
	res = DavisResistance(1e-2,0,1.5e-5)
	vᵢ = 1.0
	x, t, v, η = __solve_flat(vᵢ, track, res, V, V′)
	plot(x, v,
		xlabel = "Distance (m)",
		ylabel = "Speed (m/s",
		label = false)
end

# ╔═╡ 1e0bc6d8-d6a0-4698-9ce4-120f0d805581
begin
	x3, t3, v3, η3 = solve_flat(vᵢ, track, res, 1200)
	@show t3[end]
	plot(x3, v3)
end

# ╔═╡ d21e78b1-c193-4066-a084-5429ddba0e32
begin
	x2, t2, v2, η2 = _solve_flat(vᵢ, track, res, V)
	plot(x2, v2,
		xlabel = "Distance (m)",
		ylabel = "Speed (m/s",
		label = false)
end

# ╔═╡ 9428fe6e-7f98-4faf-bb69-b1fc8cb2462b
begin
	plot(η3, v3, label = false)
end

# ╔═╡ Cell order:
# ╠═30410070-c901-11ed-3a7b-23fe70952345
# ╟─6f33c454-b47f-4435-9623-20b4783fea11
# ╠═3b663ae4-d7ad-4e06-a3a8-88b54b52a290
# ╠═d0e4a291-ff3b-484d-959a-d185a8d71073
# ╠═425edfc2-c7d8-47c0-959f-11cd7bd58e57
# ╠═544c1c3e-cfb5-43b8-b56b-fdc7c9d30c93
# ╠═e7d708cc-5e2a-4abb-9172-8d4c045407c6
# ╟─e774edc5-4b9e-41ae-b86c-06f7b484f343
# ╟─a92c7259-91e1-4db9-840a-7f86aec18c1d
# ╠═601ae6be-7809-4c62-851f-f549dab359d1
# ╟─5d6314fa-8e7c-4e22-8afe-11603be925bc
# ╟─32dd4a90-39e7-440a-8802-bbbdd3062680
# ╟─73445f47-8029-492c-927f-2da1f9b51a75
# ╟─6a37fcbd-2ef5-49fd-8792-05c76e0720ed
# ╟─1f9e3829-f18c-433c-8b7e-928bb7f92bd3
# ╠═1e0bc6d8-d6a0-4698-9ce4-120f0d805581
# ╟─d21e78b1-c193-4066-a084-5429ddba0e32
# ╠═9428fe6e-7f98-4faf-bb69-b1fc8cb2462b
