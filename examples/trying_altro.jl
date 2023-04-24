# Code from https://roboticexplorationlab.org/Altro.jl/dev/quickstart.html
using Altro
using FiniteDiff, ForwardDiff
using LinearAlgebra
using Plots
using RobotDynamics
using StaticArrays
using TrajectoryOptimization

const TO = TrajectoryOptimization
const RD = RobotDynamics

RD.@autodiff struct Train <: RD.ContinuousDynamics end
RD.state_dim(::Train) = 2
RD.control_dim(::Train) = 1

function RD.dynamics(::Train, s, u)
    ṡ = @SVector [1 / s[2],
                  (u[1] - 1.5e-2 - 0.016e-2*s[2]^2 + 0) / s[2]]
end

RD.default_signature(::Train) = RD.StaticReturn()

RD.@autodiff struct TrainCost <: TO.CostFunction
    ρ
end
RD.state_dim(::TrainCost) = 2
RD.control_dim(::TrainCost) = 1

function RD.evaluate(cost::TrainCost, x, u)
    (abs(u[1]) + u[1]) / 2 + cost.ρ * (abs(u[1]) - u[1]) / 2
end

model = Train()
dmodel = RD.DiscretizedDynamics{RD.RK4}(model)
n, m = size(model)
N = 101
xf = 500.0
dx = xf / (N - 1)

s0 = SA_F64[0, 1.0]
sf = SA_F64[800, 1.0]

# Objective
# Q = Diagonal(SA[0.1,0.1])
# R = [1]
# Qf = Diagonal(SA[100, 100])
# obj = LQRObjective(Q, R, Qf, sf, N)
obj = TrainCost(0.0)

# Constraints
cons = ConstraintList(n, m, N)
goal = GoalConstraint(sf)
add_constraint!(cons, goal, N)

bnd = BoundConstraint(n, m, x_min = [0, 0.1])
add_constraint!(cons, bnd, 1:N-1)

ubnd = BoundConstraint(n, m, u_min = -0.25, u_max = 0.125)
add_constraint!(cons, ubnd, 1:N-1)

# Problem definition
prob = Problem(dmodel, obj, s0, xf; xf = sf, constraints = cons)

# Initialization
firstpart = fill(0.125, floor(Int, N/3))
finalpart = fill(-0.25, floor(Int, N/3))
middlepart = fill(0.05, N-1 - 2*floor(Int, N/3))
u_guess = vcat(firstpart, middlepart, finalpart)
u_guess = [[val] for val in u_guess]

initial_controls!(prob, u_guess)
rollout!(prob)

# Solving
solver = ALTROSolver(prob, iterations = 10e3)
solve!(solver)

status(solver)

S = states(solver)
U = controls(solver)

Sm = hcat(Vector.(S)...)
Um = hcat(Vector.(U)...)

plot(0:dx:xf, Sm[2,:])