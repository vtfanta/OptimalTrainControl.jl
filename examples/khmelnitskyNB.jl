### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ f9b24012-cc97-11ed-3abe-736fa4e439e9
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ 28c11473-c285-49bd-92ed-ece57bdd15f1
begin
	using RailDynamics
	using Roots
	using Plots
	using DifferentialEquations
	using DataFrames
end

# ╔═╡ f7f2edd4-357b-4fcd-9eec-3c9c1f854579
begin
trackX = [0, 16e3, 20e3, 24e3, 25e3, 28e3, 31e3, 40e3]
trackY = [0, 0, 400, 160, 160, 460, 280, 280]

mytrack = HillyTrack(trackX, trackY);
end

# ╔═╡ 1a18361a-819f-4136-9f29-01bb21e470b7
begin
	plot(mytrack)
end

# ╔═╡ 7a847255-ee9a-4336-806a-18ad5893f408
begin
	# No regeneration, i.e. ρ = 0
	vᵢ = √2
	T = 3600
	umax(v) = 0.125
	umin(v) = -0.25
	myres = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)
	# Copied from the article, have to iterate for it otherwise
	V = sqrt(63.27 * 2)
end

# ╔═╡ 2bea567f-d959-490d-b360-fd14989cb462
function getmildsegments(track, V, res, umax, umin, ρ = 0)
	midpoints = [(track.waypoints[k, :Distance] + track.waypoints[k + 1, :Distance]) / 2 for k ∈ 1:nrow(track.waypoints) - 1]
	gs = [getgradientacceleration(track, x) for x in midpoints]
	starts = [x for x in track.waypoints[1:end-1, :Distance]]
	ends = [x for x in track.waypoints[2:end, :Distance]]

	if 0 < ρ < 1
		W = find_zero(W -> ψ(res, W) - ψ(res, V) / ρ, V)
	end
	
	ret = []
	for (i, g) in enumerate(gs)
		if g ≤ resistance(res, V) ≤ umax(V) + g
			push!(ret, (starts[i], ends[i], :HoldP))
		elseif 0 < ρ < 1 &&  +g - resistance(res, W) ≥ 0
			push!(ret, (starts[i], ends[i], :HoldB))
		end
	end 
	insert!(ret, 1, (-Inf, starts[1], :HoldP))
	push!(ret, (ends[end], Inf, :HoldP))
end

# ╔═╡ a5da8052-70a6-49b9-bb5a-37aa60ae0dfe
getmildsegments(mytrack, V, myres, umax, umin)


# ╔═╡ Cell order:
# ╠═f9b24012-cc97-11ed-3abe-736fa4e439e9
# ╠═28c11473-c285-49bd-92ed-ece57bdd15f1
# ╠═f7f2edd4-357b-4fcd-9eec-3c9c1f854579
# ╠═1a18361a-819f-4136-9f29-01bb21e470b7
# ╠═7a847255-ee9a-4336-806a-18ad5893f408
# ╠═2bea567f-d959-490d-b360-fd14989cb462
# ╠═a5da8052-70a6-49b9-bb5a-37aa60ae0dfe
