# module Models
@reexport using OptimalTrainControl

export calculatecontrol!, φ, φ′, ψ, resistance, E
export getmildsegments

function calculatecontrol!(prob::TrainProblem, sol, points)
    @unpack track, umax, umin = prob
    x, v = sol.t, sol[2,:]
    u = zeros(size(x))

    pointsx = [p[1] for p in points]
    mode = nothing
    val = nothing
    for (idx,s) in enumerate(x)
        if s in pointsx
            mode = points[findfirst(s .== pointsx)][2]
        end
        if mode == :MaxP
            val = umax(v[idx])
        elseif mode == :HoldP || mode == :HoldR
            val = resistance(prob.resistance, v[idx]) - getgradientacceleration(track, s)
        elseif mode == :Coast
            val = 0
        elseif mode == :MaxB
            val = umin(v[idx])
        end
        u[idx] = val
    end
    prob.control = u
    return u
end

function getmildsegments(params::ModelParams)
    @unpack track, ρ, V, umax = params
    res = params.resistance

    if isa(track, FlatTrack)
        midpoints = [(start(track) + finish(track)) / 2]
        gs = [0]
        starts = [start(track)]
        ends = [finish(track)]
    else
        midpoints = [(track.waypoints[k, :Distance] + track.waypoints[k + 1, :Distance]) / 2 for k ∈ 1:nrow(track.waypoints) - 1]
        gs = [getgradientacceleration(track, x) for x in midpoints]
        starts = [x for x in track.waypoints[1:end-1, :Distance]]
        ends = [x for x in track.waypoints[2:end, :Distance]]
    end

	if 0 < ρ < 1
		W = find_zero(W -> ψ(res, W) - ψ(res, V) / ρ, V)
	end
	
	ret = []
	for (i, g) in enumerate(gs)
		# if g ≤ resistance(res, V) ≤ umax(V) + g
        if !(umax(V) - resistance(res, V) + g < 0) && !(-resistance(res, V) + g > 0)
			if Base.length(ret) ≥ 1 && last(ret).mode == :HoldP && starts[i] == last(ret).finish
				last(ret).finish = ends[i]
			else
				push!(ret, Segment(starts[i], ends[i], :HoldP, V))
			end
		elseif 0 < ρ < 1 &&  +g - resistance(res, W) ≥ 0
			if Base.length(ret) ≥ 1 && last(ret).mode == :HoldR && starts[i] == last(ret).finish
				last(ret).finish = ends[i]
			else
				push!(ret, Segment(starts[i], ends[i], :HoldR, W))
			end
		end
	end 
	pushfirst!(ret, Segment(-Inf, starts[1], :HoldP, V))
	push!(ret, Segment(ends[end], Inf, :HoldP, V))
end

function getmildsegments(track, V, res, umax, ρ = 0)
    @warn "This version is deprecated. Use `getmildsegments(::ModelParams)` instead."
    if isa(track, FlatTrack)
        midpoints = [(start(track) + finish(track)) / 2]
        gs = [0]
        starts = [start(track)]
        ends = [finish(track)]
    else
        midpoints = [(track.waypoints[k, :Distance] + track.waypoints[k + 1, :Distance]) / 2 for k ∈ 1:nrow(track.waypoints) - 1]
        gs = [getgradientacceleration(track, x) for x in midpoints]
        starts = [x for x in track.waypoints[1:end-1, :Distance]]
        ends = [x for x in track.waypoints[2:end, :Distance]]
    end

	if 0 < ρ < 1
		W = find_zero(W -> ψ(res, W) - ψ(res, V) / ρ, V)
	end
	
	ret = []
	for (i, g) in enumerate(gs)
		# if g ≤ resistance(res, V) ≤ umax(V) + g
        if !(umax(V) - resistance(res, V) + g < 0) && !(-resistance(res, V) + g > 0)
			if Base.length(ret) ≥ 1 && last(ret).mode == :HoldP && starts[i] == last(ret).finish
				last(ret).finish = ends[i]
			else
				push!(ret, Segment(starts[i], ends[i], :HoldP, V))
			end
		elseif 0 < ρ < 1 &&  +g - resistance(res, W) ≥ 0
			if Base.length(ret) ≥ 1 && last(ret).mode == :HoldR && starts[i] == last(ret).finish
				last(ret).finish = ends[i]
			else
				push!(ret, Segment(starts[i], ends[i], :HoldR, W))
			end
		end
	end 
	pushfirst!(ret, Segment(-Inf, starts[1], :HoldP, V))
	push!(ret, Segment(ends[end], Inf, :HoldP, V))
end

"""
Calculate the Davis formula resistant force per unit mass.

    R = a + b * v + c * v^2,

where v is the vehicle speed and a, b and c are the resistance parameters.
"""
function resistance(r::DavisResistance, u)
    if Base.length(u) > 1
        r.a + r.b * u[2] + r.c * u[2]^2
    else
        r.a + r.b * u[1] + r.c * u[1]^2
    end
end

function φ(r::DavisResistance, v)
    v * resistance(r, v)
end

function φ′(r::DavisResistance, v)
    resistance(r, v) + ψ(r, v) / v^2
end

function ψ(r::DavisResistance, u)
    if Base.length(u) > 1
        u[2]^2 * (r.b + 2r.c * u[2])
    else
        u[1]^2 * (r.b + 2r.c * u[1])
    end
end

function E(r::DavisResistance, V, v)
    ψ(r, V) / v + resistance(r, v)
end
# end # module