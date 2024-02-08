using Roots
using StaticArrays

export solve

"""
    solve(problem::TOTCProblem)

Compute `OTCSolution` of the time-optimal train control `problem`.

See also [`TOTCProblem`](@ref), [`OTCSolution`](@ref).
"""
function solve(p::TOTCProblem)
    # Time-optimal solution consists of a MaxP phase followed by a MaxB phase.
    function _odefun(s::A, p::TOTCProblem, x::T) where {T<:Real, A<:AbstractArray{T,1}}
        t, v = s

        if p.current_phase == MaxP
            u = p.train.U̅(v)
        elseif p.current_phase == MaxB
            u = p.train.U̲(v)
        end

        ds1 = 1/v
        ds2 = (u - r(p.train, v) + g(p.track, x)) / v
        SA[ds1, ds2]
    end

    if isnothing(p.track.x_gradient)
        xspan = (0., p.track.length)
    else
        xspan = (p.track.x_gradient[1], p.track.x_gradient[1] + p.track.length)
    end

    # Find maximum acceleration trajectory from start;
    # find maximum deceleration trajectory from end;
    
    s0 = SA{Float64}[0, 1.]
    maxtotcprob = TOTCProblem(p.train, p.track, MaxP)
    maxodeprob = ODEProblem(_odefun, s0, xspan, maxtotcprob)

    s0 = SA{Float64}[0, 1.]
    mintotcprob = TOTCProblem(p.train, p.track, MaxB)
    minodeprob = ODEProblem(_odefun, s0, reverse(xspan), mintotcprob)

    maxsol = OrdinaryDiffEq.solve(maxodeprob, Tsit5())
    minsol = OrdinaryDiffEq.solve(minodeprob, Tsit5())

    # find the switching point.
    x_switch = Roots.find_zero(x -> maxsol(x)[2] - minsol(x)[2], xspan)
    t_switch, v_switch = maxsol(x_switch)

    maxlastindex = searchsortedfirst(maxsol.t, x_switch)
    minfirstindex = searchsortedfirst(reverse(minsol.t), x_switch)

    tmin = reverse(minsol[1,:]) .- minimum(minsol[1,:]) .+ t_switch
    vmin = reverse(minsol[2,:])

    utot = maxsol.u[1:maxlastindex-1]
    push!(utot, [t_switch, v_switch])
    append!(utot, [[tmin[k], vmin[k]] for k=minfirstindex:length(minsol.t)])

    xtot = maxsol.t[1:maxlastindex-1]
    push!(xtot, x_switch)
    append!(xtot, reverse(minsol.t)[minfirstindex:end])

    # Construct the solution
    sol = DiffEqBase.build_solution(maxodeprob, Tsit5(), xtot, utot, retcode = ReturnCode.Success)
    x_phases = [xspan[1], x_switch]
    phases = [MaxP, MaxB]

    control = function(x)
        current_phase = phases[searchsortedlast(x_phases, x)]
        if current_phase == MaxP
            p.train.U̅(sol(x)[2])
        elseif current_phase == MaxB
            p.train.U̲(sol(x)[2])
        end
    end

    OTCSolution(sol, x_phases, phases, control)
end