using RailDynamics
using Test

@testset "RailDynamics.jl" begin
    # Write your tests here.
end

@testset "Track tests" begin
    x = [0,1,2,3,4,5]
    y = [0,1,1,2,0,1]

    dfvec = HillyTrack(x, y)

end
