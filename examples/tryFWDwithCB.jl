using Optim
using DifferentialEquations

function yeah(V)
    prob = ODEProblem(function (du, u, p, t) du[1] = u[1] end, typeof(V)[2.0], (0.0, 10.0))

    condition(u, t, int) = u[1] * one(V) - V
    affect!(int) = terminate!(int)
    cb = ContinuousCallback(condition, affect!)

    sol = solve(prob, callback = cb)
end


function f(V)
    sol = yeah(V[1])
    sum(abs2, sol[1,:])
end

Optim.minimizer(optimize(f, [10.0], BFGS(); autodiff = :forward))

