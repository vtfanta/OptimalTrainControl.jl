using Optim
using DifferentialEquations
using Plots

function f(x)
    x[1]^2 - 2x[1] + 1 + abs(x[2])
end

function g(x)
    # Hold speed and then maximal braking
    prob = ODEProblem(function (du, u, p, t)
            du[1] = inv(u[2])
            if t ≥ x[2]
                du[2] = (-0.05 - 1e-2 - 1.5e-5u[2]^2) * inv(u[2])
            else
                du[2] = 0.0
            end
        end,
        [0.0, x[1]], (0.0, 1000.0))

    condition(u, t, int) = u[2] ≤ 0.5
    affect!(int) = terminate!(int)    
    cb = DiscreteCallback(condition, affect!)

    sol = solve(prob, callback = cb)

    sum(abs2, [sol[2,end] - 1.0, sol.t[end] - 200.0])
end

o = Optim.minimizer(optimize(g, [20.0, 150.0], LBFGS(); autodiff = :forward))

prob = ODEProblem(function (du, u, p, t)
            du[1] = inv(u[2])
            if t ≥ o[2]
                du[2] = (-0.05 - 1e-2 - 1.5e-5u[2]^2) * inv(u[2])
            else
                du[2] = 0.0
            end
        end,
        [0.0, o[1]], (0.0, 1000.0))

condition(u, t, int) = u[2] ≤ 1.0
affect!(int) = terminate!(int)    
cb = DiscreteCallback(condition, affect!)

sol = solve(prob, callback = cb)

plot(sol.t, sol[2,:])