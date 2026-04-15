# use Rodas4 method to solve DAE (EETC with two states and the η state as algebraic)

using LinearAlgebra
using OptimalTrainControl
using OrdinaryDiffEq
using Plots

module myU
    export Max_u, Min_u
    import Base.max

    abstract type MaxControl end
    abstract type MinControl end

    struct Max_u{T<:Real} <: MaxControl
        P::T
        v₀::T    
    end

    (u::Max_u)(v::Real) = u.P / Base.max(u.v₀, v)

    struct Min_u{T<:Real} <: MinControl
        Q::T
        v₀::T
    end

    (u::Min_u)(v::Real) = u.Q / Base.max(u.v₀, v)
end

module MyResistance
    export DavisResistance, ψ, E

    abstract type AbstractResistance end

    struct DavisResistance{T<:Real} <: AbstractResistance
        a::T
        b::T
        c::T
    end

    (r::DavisResistance)(v::Real) = r.a + r.b * v + r.c * v^2

    ψ(r::DavisResistance, v::Real) = v^2 * (r.b + 2r.c * v)
    E(r::DavisResistance, V::Real, v::Real) = ψ(r, V) / v + r(v)
end

function rhs!(d_states, states, params, x)
    _, v, η = states
    r, umax, umin, V = params
    d_states[1] = 1/v
    d_states[2] = (umax(v) - r(v)) / v
    d_states[3] = (MyResistance.E(r, V, v) - MyResistance.E(r, V, V)) / (umax(v) - r(v)) - η
end

M = diagm([1, 1, 0])
f = ODEFunction(rhs!, mass_matrix=M)

res = MyResistance.DavisResistance(1e-2, 0.0, 1.5e-5)
u_max = myU.Max_u(1, 5)
u_min = myU.Min_u(-1, 5)
V = 25

cb = ContinuousCallback(
    (states, params, x) -> states[3],
    int -> terminate!(int),
)

dist_span = [0.0, 10530.0]
initial_states = [0.0, 1.0, 1.5]
params = (res, u_max, u_min, V)
prob = ODEProblem{true}(f, initial_states, dist_span, params)

sol = OrdinaryDiffEq.solve(prob, Rodas4(); callback=cb, dtmax=2.0)

##
plot(sol[2,:])
hline!([V])