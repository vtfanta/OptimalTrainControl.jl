using RailDynamics
using Test
using CSV
using Tables
using DataFrames

@testset "RailDynamics.jl" begin
    # Write your tests here.
end

@testset "Track tests" begin
    N = 130
    x = collect(range(0, 100, length = N))
    y = (rand(N) .- 0.5) * 2
    mat = cat(x,y, dims = 2)

    CSV.write(joinpath(pwd(), "test", "testtrack.csv"), Tables.table(mat, header=[:Distance, :Altitude]))

    htvec = HillyTrack(x, y)
    htcsv = HillyTrack(joinpath(pwd(), "test", "testtrack.csv"))
    htdfr = HillyTrack(DataFrame("Distance" => x, "Altitude" => y))

    @test all([getgrade(htvec, s) ≈ getgrade(htcsv, s) ≈ getgrade(htdfr, s)
                for s ∈ -10:0.03:110])
end
