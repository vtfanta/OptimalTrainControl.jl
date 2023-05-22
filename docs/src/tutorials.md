# Tutorial

This page shows how to use the package. 

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
```@example
using OptimalTrainControl, Plots
flattrack = FlatTrack(10e3)
prob = TrainProblem(; track = flattrack, T = 1000)
points, sol = solve!(prob)
plot(sol.t, sol[2,:]; color = modecolor(sol.t, points), lw = 2, label = false)
```

## General track

```@example
using OptimalTrainControl, Plots

# From https://doi.org/10.1109/9.867018
trackX = [0,16e3,20e3,24e3,25e3,28e3,31e3,40e3]
trackY = [0,0,400,160,160,460,280,280]/9.81
track = HillyTrack(trackX, trackY)
resistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)

T = 3600.0
u_max(v) = 0.125
u_min(v) = -0.25
vᵢ = 2.0
vf = 2.0

prob = TrainProblem(;track, resistance, T, 
    umax = u_max, umin = u_min, vᵢ, vf)

points, sol = solve!(prob)

plot(sol.t, sol[2,:]; color = modecolor(sol.t, points), lw = 2, label = false)
```