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
	using OptimalTrainControl
	using Roots
	using Plots
	using DifferentialEquations
	using DataFrames
	import ForwardDiff: derivative
	using DiffEqCallbacks
end

# ╔═╡ f7f2edd4-357b-4fcd-9eec-3c9c1f854579
begin
trackX = [0, 5e3, 6e3, 10e3]
trackY = [0, 0, 100, 100] ./ 9.81

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
		entry
		exit
	end
	function Segment(s, f, m)
		Segment(s, f, m, nothing, nothing)
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
			if length(ret) ≥ 1 && last(ret).mode == :HoldP && starts[i] == last(ret).finish
				last(ret).finish = ends[i]
			else
				push!(ret, Segment(starts[i], ends[i], :HoldP))
			end
		elseif 0 < ρ < 1 &&  +g - resistance(res, W) ≥ 0
			if length(ret) ≥ 1 && last(ret).mode == :HoldR && starts[i] == last(ret).finish
				last(ret).finish = ends[i]
			else
				push!(ret, Segment(starts[i], ends[i], :HoldB))
			end
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
	

# ╔═╡ 6e760fe2-a030-4e34-80af-c7edd35f2bc2
segs

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
	gᵦ(K) = -0.25
	aₜᵣ(ψ) = ψ ≥ 1 ? ψ - 1 : 0
	aᵦ(ψ) = ψ ≥ α ? 0 : α - ψ
	Kₛ = 63.27
	kh_track = HillyTrack(trackX, trackY * 9.81)
	P(s) = kh_track.interpolator(s)
	α = 0.0
	ψₜ = -0.3414
end

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

# ╔═╡ 1e1c6dac-163a-4c90-a264-8a6d08f48b5b
begin
	function kh_condition(out, u, t, int)
		K, ψ = u
		out[1] = K - 0.1
		# out[2] = (t ≤ segs[3].start ? 1.0 : ψ - 1.0)
		out[2] = ψ - 1.0
	end
	kh_affect!(int, idx) = terminate!(int)
	kh_cb = VectorContinuousCallback(kh_condition, nothing, kh_affect!,2)
end

# ╔═╡ 24207276-e322-4a97-a32b-6377230722c8
function kh_solve(initial_states, span)
	kh_prob = ODEProblem(khmelnitsky_odefun!, initial_states, span)
	tstops = [r[:Distance] for r in eachrow(kh_track.waypoints)]
	kh_sol = solve(kh_prob, Rosenbrock23(), callback = kh_cb, tstops = tstops,
		isoutofdomain = (u,p,t) -> u[1] ≤ 1e-2)
end

# ╔═╡ ac219b32-ee86-4f9c-87fb-9b42ca8b782f
function try_link(xstart, seg1, seg2)
		initvals = (getgradientacceleration(kh_track, seg1.finish+1) > 0) ? [Kₛ, 1.0 - 1e-3] : [Kₛ, 1.0]
		sol = kh_solve(initvals, (xstart, seg2.finish))
		@info plot(sol, idxs=(1))
		if sol.t[end] == seg2.finish
			+Inf
		elseif sol.t[end] < seg2.start
			-Inf
		else
			sol[1,end] - initvals[1]
		end
end

# ╔═╡ 6e42ba06-85cf-4edf-8d05-d94bb24d50ff
function linkage!(seg1::Segment, seg2::Segment)
	initvals = (getgradientacceleration(kh_track, seg1.finish+1) > 0) ? [Kₛ, 1.0 - 1e-3] : [Kₛ, 1.0]

	leftsol = kh_solve(initvals, (seg1.start, seg2.finish))
	rightsol = kh_solve(initvals, (seg1.start, seg2.finish))

	if leftsol.t[end] < seg2.start
		# did not even reach the second segment
		leftval = -Inf
	elseif leftsol.t[end] == seg2.finish
		# not encountering threshold value
		leftval = +Inf
	else
		leftval = leftsol[1,end] - (seg2.mode == :HoldP ? 1.0 : α)
	end

	if leftsol.t[end] == seg2.finish
		# not encountering threshold value
		rightval = +Inf
	elseif rightsol.t[end] < seg2.start
		# did not even reach the second segment
		rightval = -Inf
	else
		rightval = rightsol[1,end] - (seg2.mode == :HoldP ? 1.0 : α)
	end

	# if sign(leftval) == sign(rightval)
	# 	# no linkage possible
	# 	return leftsol, rightsol
	# end

	@show leftval, rightval

	link = find_zero(x->try_link(x,seg1,seg2), [seg1.start, seg1.finish]; xatol = 10.0)
	seg1.exit = link
end

# ╔═╡ fffd0d3b-4234-40ab-9c7a-4e9b3c470f51
link = linkage!(segs[2], segs[3])

# ╔═╡ 05b2638f-1f6d-4952-aae3-6823d1273ace
try_link(segs[2].finish,segs[2],segs[3])

# ╔═╡ 5bd25d6d-d09f-4a56-9a9e-110760922143
begin
	initstates = [Kₛ, 1.0]
	span = (4943.386367797852, segs[3].finish)
	kh_sol = kh_solve(initstates, span)
	solK = kh_sol[1,:]
	Ψ = kh_sol[2,:]
	@show kh_sol.t[end]
end

# ╔═╡ 1acab50f-91fe-4c8a-a7be-c5a34b45a81f
begin
	plot(kh_track)
	plot!(twinx(), kh_sol.t, solK)
end

# ╔═╡ b4194b03-d200-4395-bbc7-936b475983b2
plot(kh_sol.t, Ψ)

# ╔═╡ e5065f63-49ad-45ac-8161-c9da1091a950
kh_sol.t[end]

# ╔═╡ a2bf3e18-a114-4198-8005-478cff581007
kh_sol.retcode

# ╔═╡ Cell order:
# ╟─7ad7cae2-dee3-463b-a5c7-8754b3f15d19
# ╟─f9b24012-cc97-11ed-3abe-736fa4e439e9
# ╟─28c11473-c285-49bd-92ed-ece57bdd15f1
# ╠═f7f2edd4-357b-4fcd-9eec-3c9c1f854579
# ╟─1a18361a-819f-4136-9f29-01bb21e470b7
# ╠═7a847255-ee9a-4336-806a-18ad5893f408
# ╟─2bea567f-d959-490d-b360-fd14989cb462
# ╠═a5da8052-70a6-49b9-bb5a-37aa60ae0dfe
# ╠═6e760fe2-a030-4e34-80af-c7edd35f2bc2
# ╟─81a21399-759e-46d8-bfc7-bbd1105d69c8
# ╟─fedce6c9-d6ad-45df-9fbc-917a7df55122
# ╟─b05d5976-32f0-44a8-acc5-89d7cb06d00f
# ╟─a1ef92ac-f6d9-4620-83d8-95e55935d2e7
# ╠═ac219b32-ee86-4f9c-87fb-9b42ca8b782f
# ╟─6e42ba06-85cf-4edf-8d05-d94bb24d50ff
# ╠═fffd0d3b-4234-40ab-9c7a-4e9b3c470f51
# ╠═7620ecaa-6a14-46c1-972e-b9c27c480b1f
# ╠═1eec1a0b-7f55-4891-82a3-05a75e625f0b
# ╠═ff997eca-3ce9-4dc0-9c2e-f1d6fc34b798
# ╟─1e1c6dac-163a-4c90-a264-8a6d08f48b5b
# ╠═24207276-e322-4a97-a32b-6377230722c8
# ╠═05b2638f-1f6d-4952-aae3-6823d1273ace
# ╠═5bd25d6d-d09f-4a56-9a9e-110760922143
# ╠═1acab50f-91fe-4c8a-a7be-c5a34b45a81f
# ╠═b4194b03-d200-4395-bbc7-936b475983b2
# ╠═e5065f63-49ad-45ac-8161-c9da1091a950
# ╠═a2bf3e18-a114-4198-8005-478cff581007
