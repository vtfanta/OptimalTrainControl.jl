# OptimalTrainControl.jl

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://vtfanta.github.io/OptimalTrainControl.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://vtfanta.github.io/OptimalTrainControl.jl/dev/)
[![Build Status](https://github.com/vtfanta/OptimalTrainControl.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vtfanta/OptimalTrainControl.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/vtfanta/OptimalTrainControl.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/vtfanta/OptimalTrainControl.jl)

Package for energy-efficient speed profile planning for a single train.

# Installation

Currently, the package can be installed by entering
```julia
using Pkg; Pkg.add(url="https://github.com/vtfanta/OptimalTrainControl.jl")
```
to the Julia REPL.

# Documentation

The package reference and documentation is available [here](https://vtfanta.github.io/OptimalTrainControl.jl/dev/)
with a selection of presented examples. A most basic example is also presented below.

This package was developed as a part of master's thesis which can be read [here](https://dspace.cvut.cz/handle/10467/109765?locale-attribute=enhttps://dspace.cvut.cz/handle/10467/109765?locale-attribute=en) (along with assessments) or as a [PDF](https://dspace.cvut.cz/bitstream/handle/10467/109765/F3-DP-2023-Fanta-Vit-DP_fantavit.pdf?sequence=-1&isAllowed=y).

# Basic Example
Consider the problem of finding the optimal speed profile on a flat track of length $10\ \mathrm{km}$. The total journey time is
$1000\ \mathrm{s}$ and the other parameters are set to their default values.
The parameter values are passed as keyword arguments to the `TrainProblem` struct:
```julia
using OptimalTrainControl, Plots 
flattrack = FlatTrack(10e3)
prob = TrainProblem(; track = flattrack, T = 1000)
```
We can now solve the problem and plot the results:
```julia
points, sol = solve!(prob)
plot(sol.t, sol[2,:]; color = modecolor(sol.t, points), lw = 2, label = false,
    xlabel = "Distance (m)", ylabel = "Speed (m/s)") 
```
![Flat track solution](/examples/flattrack_solution.png)


