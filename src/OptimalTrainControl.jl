module OptimalTrainControl

# Type definitions and constructors
include("types.jl")

# Track-related functions
include("track.jl")

# Solving time-optimal (minimal time) train control
include("time_optimal.jl")

# Solving energy-optimal train control on a flat track
include("optimal_flat.jl")

# Simulation-related functions for regular phases and η calculation
include("simulation.jl")

# Solving energy-optimal train control on a track with gradients
include("energy_optimal.jl")

# Plotting-related code
include("plot_recipes.jl")

# TTOBench track loading
include("ttobench_loading.jl")

# Precompilation
# include("precompile.jl")

end
