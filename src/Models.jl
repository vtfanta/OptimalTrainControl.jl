# module Models
@reexport using OptimalTrainControl

export albrecht_odefun!
export play, controllaw, calculatecontrol!, φ, φ′, ψ, resistance, E
export getmildsegments

"""
Give the differential equations given in Albrecht et al. 2016, equations (1) and (2).
First state is time (second), second state is speed (metre per second).
"""
function albrecht_odefun!(du, u, p, t)
    du[1] = inv(u[2])
    du[2] = (p.u(u, p, t) - p.r(u, p, t) + p.g(u, p, t)) * inv(u[2])
end

function getmildsegments(params::NewModelParams)
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
    @warn "This version is deprecated. Use `getmildsegments(::NewModelParams)` instead."
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



function play(s::BasicScenario, initialvalues)
    X = 100.0
    prob = ODEProblem(odefun(s), initialvalues, (0.0, X))
    solve(prob)
end

function play(s::MinimalTimeScenario)
    if isnothing(s.controllaw)
        error("The scenario does not have a control law. Try running `calculatecontrols!(scenario)`.")
    end

    condition(u, t, int) = u[1] ≤ s.finalvalues[1]
    affect!(int) = terminate!(int)    
    cb = DiscreteCallback(condition, affect!)

    prob = ODEProblem(odefun(s), s.initialvalues, (0.0, length(s.track)))
    sol = solve(prob, 
                dense = false, 
                callback = cb,
                dtmax = 200)
end

function odefun(s::OptimalScenario)
    m = s.model
    return function (du, u, p, t)
        v = u[2]
        μ₂ = u[3]

        η = μ₂ * inv(v) - 1
        ζ = η + 1 - m.ρ

        μ₁ = -ψ(m.resistance, s.V)

        du[1] = inv(v)
        du[2] = (m.controllaw(t, v, η, ζ) - resistance(m.resistance, v) + 
            getgradientacceleration(s, t)) * inv(v)
        du[3] = μ₁ / v^2 + μ₂ * du[2] * inv(v) + 
            μ₂ * ψ(m.resistance, v) * inv(v)^3 - 
            (η > 0 ? η : 0) * derivative(m.maxcontrol, v) + 
            (ζ < 0 ? -ζ : 0) * derivative(m.mincontrol, v)
    end
end

function odefun(s::MinimalTimeScenario)
    m = s.model
    return function (du, u, p, t)
        du[1] = (s.controllaw(u, t) - resistance(m.resistance, u[1]) + getgradientacceleration(s, t)) * inv(u[1])
    end
end

function odefun(s::BasicScenario)
    m = s.model
    return function (du, u, p, t)
        du[1] = inv(u[2])
        du[2] = (m.maxcontrol(u[2]) - resistance(m.resistance, u) + getgradientacceleration(s, t)) * inv(u[2])
    end
end

function calculatecontrol!(s::OptimalScenario)
    # TODO
end

function calculatecontrol!(s::MinimalTimeScenario)
    function _maxthrottle!(du, u, p, t)
        du[1] = (s.model.maxcontrol(u[1]) - resistance(s.model.resistance, u) + 
            getgradientacceleration(s, t)) * inv(u[1])
    end
    function _maxbrake!(du, u, p, t)
        du[1] = (s.model.mincontrol(u[1]) - resistance(s.model.resistance, u) + 
            getgradientacceleration(s, t)) * inv(u[1])
    end

    @info "Calculating the time-optimal strategy"

    throttleprob = ODEProblem(_maxthrottle!, s.initialvalues, (0.0, length(s.track)))
    brakeprob = ODEProblem(_maxbrake!, s.finalvalues, (length(s.track), 0))

    throttlesol = solve(throttleprob)
    brakesol = solve(brakeprob)

    throttleinterpolator = 
        CubicSplineInterpolator(throttlesol.t, throttlesol[1,:], WeakBoundaries())
    brakeinterpolator = 
        CubicSplineInterpolator(reverse(brakesol.t), reverse(brakesol[1,:]), WeakBoundaries())

    switchingpoint = 
        find_zero(x -> throttleinterpolator(x) - brakeinterpolator(x), 1.0)

    # @show switchingpoint

    s.controllaw = 
        function (u, t) 
            # @show t
            if t ≤ switchingpoint
                s.model.maxcontrol(u[1])
            else
                s.model.mincontrol(u[1])
            end
        end
    throttlesol, brakesol, throttleinterpolator, brakeinterpolator
end

function controllaw(m::Model, u, p, t)
    # All gas, no brakes
    m.maxcontrol
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