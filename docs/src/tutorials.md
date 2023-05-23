# Tutorial

This page shows how to use the OptimalTrainControl.jl. 

## Track creation
A flat track instance can be created very simply since it is defined only by its length. A flat
track of a length $10\ \mathrm{km}$ can be created with
```julia
using OptimalTrainControl
flattrack = FlatTrack(10e3)
```

A general track with a nontrivial terrain is defined with the use of two vectors which represent
the distances and height values at the corresponding distance (in metres):
```@example
using OptimalTrainControl, Plots
trackX = [0, 2, 3, 5, 7]
trackY = [0, 0.2, -0.1, -0.1, 0.1]
track = HillyTrack(trackX, trackY)
plot(track)
```

## Flat track
Consider now the problem of finding the optimal speed profile on a flat track defined above. The total journey time is
$1000\ \mathrm{s}$ and the other parameters are set to their default values. The default values are
`vᵢ = 1.0`, `vf = 1.0`, `umax = v -> 3 / max(5, v)`, `umin = v -> -3 / max(5, v)`, `ρ = 0`, ``r(v)=0.01 + 1.5\cdot 10^{-5}v^2``.
The parameter values are passed as keyword arguments to the `TrainProblem` struct:
```julia
prob = TrainProblem(; track = flattrack, T = 1000)
```
We can now solve the problem and plot the results:
```@example
using OptimalTrainControl, Plots # hide
flattrack = FlatTrack(10e3) # hide
prob = TrainProblem(; track = flattrack, T = 1000) # hide
points, sol = solve!(prob)
plot(sol.t, sol[2,:]; color = modecolor(sol.t, points), lw = 2, label = false,
    xlabel = "Distance (m)", ylabel = "Speed (m/s)") 
```

## General track
The previous example involved a flat track. However, the Earth, as some of us still believe, is not as flat as it would seem. The following example shows calculation of optimal speed profile for the case of a track containing steep sections.

The track gradient profile is taken from [this article](https://doi.org/10.1109/9.867018) written by Eugene Khmelnitsky. Let's visualise it:
```@example
using OptimalTrainControl, Plots # hide
trackX = [0,16e3,20e3,24e3,25e3,28e3,31e3,40e3]
trackY = [0,0,400,160,160,460,280,280]/9.81
track = HillyTrack(trackX, trackY)
plot(track)
```
The steep sections in the middle of the track are what makes this track a bit trickier to
plan around. The problem data are 
```math
U_+(v) = 0.125,\quad U_-(v) = -0.25,\\
r(v) = 1.5\cdot10^{-2} + \frac{0.127\cdot10^{-2}}{\sqrt2}v+\frac{0.016\cdot10^{-2}}{2}v^2,\\
v_i = v_f = 2,\\
T = 3600.
```
Putting it all together:
```@example
using OptimalTrainControl, Plots # hide
trackX = [0,16e3,20e3,24e3,25e3,28e3,31e3,40e3] # hide
trackY = [0,0,400,160,160,460,280,280]/9.81 # hide
track = HillyTrack(trackX, trackY) # hide
resistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)
T = 3600.0
u_max(v) = 0.125
u_min(v) = -0.25
vᵢ = 2.0
vf = 2.0

prob = TrainProblem(;track, resistance, T, 
    umax = u_max, umin = u_min, vᵢ, vf)

points, sol = solve!(prob)

plot(sol.t, sol[2,:]; color = modecolor(sol.t, points), lw = 2, label = false,
    xlabel = "Distance (m)", ylabel = "Speed (m/s)")
plot!(twinx(), track; xlabel = "")
```

Apart from the `sol` structure (which is taken directly from the output of `DifferentialEquaionts.jl`) the `points`
list contains the switching locations between individual control modes. The nomenclature
is taken from [this comprehensive summary](https://doi.org/10.1016/j.trb.2015.07.023) by A. Albrecht et al.
The used control mode abbreviations are: `MaxP` for maximum acceleration, `MaxB` for maximum braking,
`HoldP` for holding constant speed using positive control, `HoldR` for holding constant speed using
negative control, `Coast` for coasting (no braking or accelerating).
```@repl
using OptimalTrainControl, Plots # hide
trackX = [0,16e3,20e3,24e3,25e3,28e3,31e3,40e3] # hide
trackY = [0,0,400,160,160,460,280,280]/9.81 # hide
track = HillyTrack(trackX, trackY) # hide
resistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2) # hide
T = 3600.0 # hide
u_max(v) = 0.125 # hide
u_min(v) = -0.25 # hide
vᵢ = 2.0 # hide
vf = 2.0 # hide
prob = TrainProblem(;track, resistance, T, umax = u_max, umin = u_min, vᵢ, vf) # hide
points, sol = solve!(prob) # hide
points
```

## Recalculation
Suppose that circumstances have caused the train to deviate from the precomputed
optimal speed profile. The optimisation can be also performed on a `subtrack` of the
original one with different boundary conditions:
```@example
using OptimalTrainControl, Plots # hide
trackX = [0,16e3,20e3,24e3,25e3,28e3,31e3,40e3] # hide
trackY = [0,0,400,160,160,460,280,280]/9.81 # hide
track = HillyTrack(trackX, trackY) # hide
resistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2) # hide
T = 2600.0
u_max(v) = 0.125 # hide
u_min(v) = -0.25 # hide
vᵢ = 12.0
vf = 2.0
prob = TrainProblem(;track = subtrack(track, 11e3), resistance, T, umax = u_max, umin = u_min, vᵢ, vf)

points, sol = solve!(prob)

plot(sol.t, sol[2,:]; color = modecolor(sol.t, points), lw = 2, label = false,
    xlabel = "Distance (m)", ylabel = "Speed (m/s)")
plot!(twinx(), subtrack(track, 11e3); xlabel = "")
```