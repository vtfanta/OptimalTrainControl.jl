### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ f9b24012-cc97-11ed-3abe-736fa4e439e9
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ 7ad7cae2-dee3-463b-a5c7-8754b3f15d19
Pkg.add("DiffEqCallbacks")

# ╔═╡ 28c11473-c285-49bd-92ed-ece57bdd15f1
begin
	using RailDynamics
	using Roots
	using Plots
	using DifferentialEquations
	using DataFrames
	import ForwardDiff: derivative
	using DiffEqCallbacks
end

# ╔═╡ f7f2edd4-357b-4fcd-9eec-3c9c1f854579
begin
trackX = [0, 16e3, 20e3, 24e3, 25e3, 28e3, 31e3, 40e3]
trackY = [0, 0, 400, 400, 400, 400, 400, 400] ./ 9.81

mytrack = HillyTrack(trackX, trackY);
end

# ╔═╡ 1a18361a-819f-4136-9f29-01bb21e470b7
begin
	plot(mytrack)
end

# ╔═╡ 7a847255-ee9a-4336-806a-18ad5893f408
begin
	# No regeneration, i.e. ρ = 0
	T = 3600
	umax(v) = 0.125
	umin(v) = -0.25
	ρ = 0
	myres = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)
	# Copied from the article, have to iterate for it otherwise
	V = sqrt(63.27 * 2)
end

# ╔═╡ 81a21399-759e-46d8-bfc7-bbd1105d69c8
begin
	mutable struct Segment
		start
		finish
		mode
		s
		v
		η
	end
	function Segment(s, f, m)
		Segment(s, f, m, nothing, nothing, nothing)
	end
	Base.show(io::IO, s::Segment) = show(io,"[$(s.start)..$(s.mode)..$(s.finish)]")
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
			push!(ret, Segment(starts[i], ends[i], :HoldP))
		elseif 0 < ρ < 1 &&  +g - resistance(res, W) ≥ 0
			push!(ret, Segment(starts[i], ends[i], :HoldB))
		end
	end 
	insert!(ret, 1, Segment(-Inf, starts[1], :HoldP))
	push!(ret, Segment(ends[end], Inf, :HoldP))
end

# ╔═╡ a5da8052-70a6-49b9-bb5a-37aa60ae0dfe
begin
	segs = getmildsegments(mytrack, V, myres, umax, umin)
	S2, S3 = segs[2], segs[3]
end
	


# ╔═╡ fedce6c9-d6ad-45df-9fbc-917a7df55122
function mycontrol(u, p, x)
	η = u[3] / u[2] - 1.0
	# if η ≥ 0
		umax(u[2])
	# elseif -1 ≤ η ≤ 0
	# 	0
	# else
	# 	umin(u[2])
	# end
end

# ╔═╡ b05d5976-32f0-44a8-acc5-89d7cb06d00f
function odefun!(du, u, p, x)
        t, v, μ₂ = u
    
        η = μ₂ / v - 1.0
        ζ = η + 1.0 - ρ
    
        π₁ = η > 0. ? η : 0.
        π₂ = ζ ≥ 0. ? 0. : -ζ

        # @show v, μ₂, x, mycontrol(v, μ₂, x)
    
        du[1] = inv(v)
        du[2] = (mycontrol(u, p, x) - resistance(myres, v) + 
            getgradientacceleration(mytrack, x)) * inv(v)
        du[3] = -ψ(myres, V) / v^2 + μ₂ * du[2] / v + 
            μ₂ * derivative(x -> resistance(myres, x), v) / v - 
            π₁ * derivative(x -> umax(x), v) 
            π₂ * derivative(x -> umin(x), v)
end

# ╔═╡ a1ef92ac-f6d9-4620-83d8-95e55935d2e7
begin
	sspan = (15750, S3.finish)
	vᵢ = 0.1
	μ₂ᵢ = vᵢ * (E(myres, V, vᵢ) - E(myres, V, V)) / (umax(vᵢ) - resistance(myres, vᵢ) + getgradientacceleration(mytrack, sspan[1])) + vᵢ
	prob = ODEProblem(odefun!, [0.0, V, V], sspan)

	function condition(out, u, t, p)
		out[1] = t < S3.start ? 1 : abs(u[2] - u[3]) - 1e-3
		out[2] = u[2] - 1e-2
	end
	affect!(int, idx) = terminate!(int)
	saved_values = SavedValues(Float64, Float64)
	cb1 = SavingCallback((u,t,integrator)->mycontrol(u,integrator.p,t), saved_values)
	cb2 = VectorContinuousCallback(condition, affect!, 2)
	cbs = CallbackSet(cb1, cb2)
	sol = solve(prob, Rodas4P(); alg_hints = [:stiff], callback = cbs)

	v = sol[2,:]
	s = sol.t
	η = sol[3,:] ./ sol[2,:] .- 1
end

# ╔═╡ 1eec1a0b-7f55-4891-82a3-05a75e625f0b
begin
	w(K) = 1.5e-2 + 0.127e-2√abs(K) + 0.016e-2 * K
	R(K) = abs(K)^(3/2) * derivative(w, K)
	gₜᵣ(K) = 0.125
	gᵦ(K) = 0.25
	aₜᵣ(ψ) = ψ ≥ 1 ? ψ - 1 : 0
	aᵦ(ψ) = ψ ≥ α ? 0 : α - ψ
	Kₛ = 63.27
	kh_track = HillyTrack(trackX, trackY * 9.81)
	P(s) = kh_track.interpolator(s)
	α = 0.0
	ψₜ = -0.3414
end

# ╔═╡ 1ec8aaea-d0f7-4380-824c-d4d29a4dec76
plot(x->derivative(P,x), start(kh_track), finish(kh_track))

# ╔═╡ 7620ecaa-6a14-46c1-972e-b9c27c480b1f
function khmelnitsky_u(ψ, K)
	if ψ ≥ 1.0
		umax(K)
	elseif α ≤ ψ ≤ 1.0
		0
	else
		umin(K)
	end
end

# ╔═╡ ff997eca-3ce9-4dc0-9c2e-f1d6fc34b798
function khmelnitsky_odefun!(du, u, p, t)
	K, ψ = u
	
	du[1] = khmelnitsky_u(ψ, K) - w(K) - derivative(P, t)
	# du[1] = gₜᵣ(K) - w(K) - derivative(P, t)
	du[2] = (ψ * R(K) - R(Kₛ)) / abs(K)^(3/2) - aₜᵣ(ψ) * derivative(gₜᵣ, K)
		- aᵦ(ψ) * derivative(gᵦ, K)
	# du[2] = ψₜ / (2abs(K))^(3/2) + ψ * derivative(w, K) - aₜᵣ(ψ) * derivative(gₜᵣ, K) - aᵦ(ψ) * derivative(gᵦ, K)
	# @show du[2]
end

# ╔═╡ 6e42ba06-85cf-4edf-8d05-d94bb24d50ff
function linkage(seg1::Segment, seg2::Segment)
	function condition!(out, u, t, int)
		K, Ψ = u
		out[1] = K - 0.1
		out[2] = t < seg2.start ? 1.0 : Ψ - (seg2.mode == :HoldP ? 1.0 : α)
	end
	
	affect!(int, idx) = terminate!(int)
	cb = VectorContinuousCallback(condition!, affect!, 2)
	
	function try_link(xstart)
		prob = ODEProblem(khmelnitsky_odefun!, initvals, (xstart, seg2.finish);
		callback = cb)
		sol = solve(prob, Rosenbrock23())

		if sol.t[end] == seg2.finish
			+Inf
		elseif sol.t[end] < seg2.start
			-Inf
		else
			sol[1,end] - initvals[1]
		end
	end

	if seg1.finish == seg2.start
		# no linkage possible on adjacent segments
		return false
	end

	initvals = [Kₛ, 1.0]
	leftprob = ODEProblem(khmelnitsky_odefun!, initvals, (seg1.start, seg2.finish);
		callback = cb)
	rightprob = ODEProblem(khmelnitsky_odefun!, initvals, (seg1.finish, seg2.finish);
		callback = cb)

	leftsol = solve(leftprob, Rosenbrock23())
	rightsol = solve(rightprob, Rosenbrock23())

	if leftsol.t[end] < seg2.start
		# did not even reach the second segment
		leftval = -Inf
	elseif leftsol.retcode == :Success
		# not encountering threshold value
		leftval = +Inf
	end

	if leftsol.retcode == :Success
		# not encountering threshold value
		rightsol = +Inf
	end

	if leftsol == rightsol
		# no linkage possible
		return false
	end

	link = find_zero(try_link, [seg1.start, seg1.finish])
	
end

# ╔═╡ fffd0d3b-4234-40ab-9c7a-4e9b3c470f51
linkage(S2, S3)

# ╔═╡ 1e1c6dac-163a-4c90-a264-8a6d08f48b5b
begin
	function kh_condition(out, u, t, int)
		K, ψ = u
		out[1] = K - 0.1
		out[2] = t < S3.start ? 1.0 : ψ - 1
	end
	kh_affect!(int, idx) = terminate!(int)

	U = SavedValues(Float64, Float64)
	saving_cb = SavingCallback((u,t,int) -> khmelnitsky_u(u[2],u[1]), U)
	kh_cb1 = VectorContinuousCallback(kh_condition, kh_affect!, 2)
	kh_cb = CallbackSet(saving_cb, kh_cb1)
end

# ╔═╡ 24207276-e322-4a97-a32b-6377230722c8
begin
	kh_span = (15750, S3.finish)
	initial_states = [Kₛ, 1.0]
	kh_prob = ODEProblem(khmelnitsky_odefun!, initial_states, kh_span)
	U;
	kh_sol = solve(kh_prob, Rosenbrock23(), callback = kh_cb,
		isoutofdomain = (u,p,t) -> u[1] ≤ 1e-2)

	solK, Ψ = kh_sol[1,:], kh_sol[2,:]
end

# ╔═╡ 1acab50f-91fe-4c8a-a7be-c5a34b45a81f
begin
	plot(kh_track)
	plot!(twinx(), kh_sol.t, solK)
end

# ╔═╡ da815f6f-fe07-4808-b82a-d85fc02f03e4
plot(kh_sol.t, Ψ)

# ╔═╡ b4194b03-d200-4395-bbc7-936b475983b2
kh_sol.t[end]

# ╔═╡ Cell order:
# ╠═7ad7cae2-dee3-463b-a5c7-8754b3f15d19
# ╠═f9b24012-cc97-11ed-3abe-736fa4e439e9
# ╠═28c11473-c285-49bd-92ed-ece57bdd15f1
# ╠═f7f2edd4-357b-4fcd-9eec-3c9c1f854579
# ╠═1a18361a-819f-4136-9f29-01bb21e470b7
# ╠═7a847255-ee9a-4336-806a-18ad5893f408
# ╟─2bea567f-d959-490d-b360-fd14989cb462
# ╠═a5da8052-70a6-49b9-bb5a-37aa60ae0dfe
# ╟─81a21399-759e-46d8-bfc7-bbd1105d69c8
# ╟─fedce6c9-d6ad-45df-9fbc-917a7df55122
# ╟─b05d5976-32f0-44a8-acc5-89d7cb06d00f
# ╟─a1ef92ac-f6d9-4620-83d8-95e55935d2e7
# ╠═1ec8aaea-d0f7-4380-824c-d4d29a4dec76
# ╠═6e42ba06-85cf-4edf-8d05-d94bb24d50ff
# ╠═fffd0d3b-4234-40ab-9c7a-4e9b3c470f51
# ╟─7620ecaa-6a14-46c1-972e-b9c27c480b1f
# ╠═1eec1a0b-7f55-4891-82a3-05a75e625f0b
# ╠═ff997eca-3ce9-4dc0-9c2e-f1d6fc34b798
# ╠═1e1c6dac-163a-4c90-a264-8a6d08f48b5b
# ╠═24207276-e322-4a97-a32b-6377230722c8
# ╠═1acab50f-91fe-4c8a-a7be-c5a34b45a81f
# ╠═da815f6f-fe07-4808-b82a-d85fc02f03e4
# ╠═b4194b03-d200-4395-bbc7-936b475983b2
