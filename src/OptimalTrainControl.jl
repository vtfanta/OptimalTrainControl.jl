module OptimalTrainControl

# Type definitions and constructors
include("types.jl")

# Track-related functions
include("track.jl")

# Solving time-optimal (minimal time) train control
include("time_optimal.jl")

# Solving energy-optimal train control on a flat track
include("optimal_flat.jl")

# Plotting-related code
include("plot_recipes.jl")

# TTOBench track loading
include("ttobench_loading.jl")

# Precompilation
include("precompile.jl")

end
