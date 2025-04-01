using OptimalTrainControl

#TODO can V be put inside the ODEProblem parameters?
"""
Create condition for continuous callback ODE solver function (root finding) for the transition
between MaxP and Coast phases. The factory is needed because the calculation of η depends on the 
speed parameter V.
"""
function make_cond_power2coast(p::EETCProblem, V)
    return function cond_power2coast(s, x, int)
        η = calculate_η(s, x, int.p, V)
        η
    end
end

"""
Create condition for continuous callback ODE solver function (root finding) for the transition
between Coast and Brake phases. The factory is needed because the calculation of η depends on the 
speed parameter V.
"""
function make_cond_coast2brake(p::EETCProblem, V)
    return function cond_coast2brake(s, x, int)
        η = calculate_η(s, x, int.p, V)
        η - (int.p.train.ρ - 1.)
    end
end

"""
Create affect! for continuous callback ODE solver function (root finding) for the transition
to the Coast phase. The factory is needed because calculation of E depends on the
speed parameter V.

Augment the vector of the E coefficients which depend on current value of η to assure
continuity.
"""
function make_aff_2coast!(p::EETCProblem, V)
    return function aff_2coast!(int)
        @info "To Coast!"
        if int.p.current_phase == MaxP
            int.p.current_phase = Coast
        elseif int.p.current_phase == MaxB
            int.p.current_phase = Coast
            # add the last E to the E array (since now there are Fs because we are in MaxB)
            push!(int.p.Es, -E(int.p.train, V, int.u[2]))
        end
    end
end

"""
Create affect! for continuous callback ODE solver function (root finding) for the transition
to the Brake phase. The factory is needed because the calculation of E depends on the
speed parameter W which itself depends on V, train resistance and ρ.

Augment the vector of the E coefficients which depend on current value of η to assure
continuity.
"""
function make_aff_2brake!(p::EETCProblem, V)
    W = p.train.ρ > 0. ? Roots.find_zero(v -> -ψ(p.train, V) + p.train.ρ*ψ(p.train, v), V) : V
    function aff_2brake!(int)
        @info "to Brake!"
        int.p.current_phase == Coast ? int.p.current_phase = MaxB : nothing

        # add the last F to the E array
        # ζ = 0 -> F = -ρ * E(W, v)
        push!(int.p.Es, -int.p.train.ρ * E(int.p.train, W, int.u[2]))
    end
end

"""
Create affect! for continuous callback ODE solver function (root finding) for the transition
to the Power phase. The factory is not necessarily needed, but other affect! functions are factories.
"""
aff_2power!(int) = int.p.current_phase == Coast ? int.p.current_phase = MaxP : nothing

"""
Right-hand side of the ODE system for the EETC problem. The function calculates the derivatives
and returns them in a StaticArray. The function assumes presence in a regular phase (MaxP, Coast, MaxB), so
it is used only to link singular phases.
"""
function _odefun(s::A, p::EETCProblem, x::T) where {T<:Real, A<:AbstractArray{T,1}}
    t, v = s

    if p.current_phase == MaxP
        u = p.train.U̅(v)
    elseif p.current_phase == Coast
        u = 0.
    elseif p.current_phase == MaxB
        u = p.train.U̲(v)
    end
    # @show u

    ds1 = 1/v
    ds2 = (u - OptimalTrainControl.r(p.train, v) + OptimalTrainControl.g(p.track, x)) / v
    SA[ds1, ds2]
end


# saving η (ζ shifted)
function calculate_η(s, x, p::EETCProblem, V)
    v = s[2]
    
    if p.current_phase == MaxP
        val = (E(p.train, V, v) + last(p.Es)) / 
            (p.train.U̅(v) - OptimalTrainControl.r(p.train, v) + g(p.track, x))
        # @show val, x, p.current_phase
    elseif p.current_phase == Coast
        val = (E(p.train, V, v) + last(p.Es)) / 
            (- OptimalTrainControl.r(p.train, v) + g(p.track, x))
    elseif p.current_phase == MaxB
        # have to calculate ζ, the Fs are also in the EETC.Es array
        val = (p.train.ρ*E(p.train, W,v) + last(p.Es)) /
            (p.train.U̲(v) - OptimalTrainControl.r(p.train, v) + g(p.track, x))
        # η = ζ + ρ - 1
        val += p.train.ρ - 1.
    end
    return val
end

# calculate E coefficients for individual gradient transitions
function make_update_Es(V)
    return function update_Es!(s, x, int)
        p = int.p
        η_current = calculate_η(s, x-1e-4, p, V) # just before gradient change
        # @show η_current, p.current_phase
        g_prev = g(p.track, x-0.01)
        g_next = g(p.track, x+0.01)
        if p.current_phase == MaxP || p.current_phase == Coast  # save E
            push!(p.Es, (g_next - g_prev) * η_current + last(p.Es))
            # @show calculate_η(s, x+0.1, p, V) 
        elseif p.current_phase == MaxB  # save F
            push!(p.Es, (g_next - g_prev) * (η_current + 1 - p.train.ρ) + last(p.Es))
        end
    end
end