using RailDynamics
using Test
using CSV
using Tables
using DataFrames

@testset "Track methods" begin
    N = 130
    x = collect(range(0, 100, length = N))
    y = (rand(N) .- 0.5) * 2
    mat = cat(x,y, dims = 2)

    CSV.write(joinpath(pwd(), "testtrack.csv"), Tables.table(mat, header=[:Distance, :Altitude]))

    htvec = HillyTrack(x, y)
    htcsv = HillyTrack(joinpath(pwd(), "testtrack.csv"))
    htdfr = HillyTrack(DataFrame("Distance" => x, "Altitude" => y))

    @test all([getgrade(htvec, s) ≈ getgrade(htcsv, s) ≈ getgrade(htdfr, s)
                for s ∈ -10:0.03:110])

    flattrack = FlatTrack(100.)
    @test start(flattrack) ≈ 0
    @test finish(flattrack) ≈ 100
    @test getgrade(flattrack, (finish(flattrack) - start(flattrack)) / 2) ≈ 0
end

@testset "Segments, ρ = 0" begin
    flattrack = FlatTrack(2e3)
    V = 10
    res = DavisResistance(1e-2, 0, 1.5e-5)
    u_max(v) = 1/max(5, v)
    segs = getmildsegments(flattrack, V, res, u_max)
    @test length(segs) == 3 # Flat track only has 3 segments (initial, the track, final)
    @test segs[2].start ≈ start(flattrack)
    @test segs[2].finish ≈ finish(flattrack)

    # Example for the https://doi.org/10.1016/j.automatica.2009.07.028
    trackX = [0, 5e3, 5.6e3, 6e3, 6.5e3, 6.8e3, 10e3]
    γs = [-0.075, -0.22, -0.27, -0.15, -0.2, -0.09]
    res = DavisResistance(6.75e-3, 0, 5e-5)
    V = 20
    u_max(v) = 3 / v
    hillytrack = HillyTrack(trackX, γs)
    segs = getmildsegments(hillytrack, V, res, u_max)
    @test length(segs) == 4
    @test segs[2].finish ≈ 5e3
    @test segs[3].start ≈ 6.8e3
end

@testset "Segments, 0 < ρ < 1" begin
    # From https://doi.org/10.1109/9.867018
    trackX = [0,16e3,20e3,24e3,25e3,28e3,31e3,40e3]
    trackY = [0,0,400,160,160,460,280,280]/9.81
    track = HillyTrack(trackX, trackY)
    myresistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)
    V = sqrt(2 * 65.43)
    ρ = 0.7

    segs = getmildsegments(track, V, myresistance, v -> 0.125, ρ)
    @test length(segs) == 7
    @test segs[3].mode == :HoldR && segs[5].mode == :HoldR
end

@testset "FlatTrack" begin
    track = FlatTrack(10e3)
    T = 1600
    prob = TrainProblem(;track, T)
    atol = 5

    chain, sol = solve!(prob; atol)
    @test abs(sol[1,end] - T) ≤ atol
    @test chain[1][1] == start(track)
    @test chain[end][1] == finish(track)
end

@testset "Linking, ρ = 0" begin
    # From https://doi.org/10.1109/9.867018
    trackX = [0,16e3,20e3,24e3,25e3,28e3,31e3,40e3]
    trackY = [0,0,400,160,160,460,280,280]/9.81
    track = HillyTrack(trackX, trackY)
    myresistance = DavisResistance(1.5e-2, 0.127e-2/sqrt(2), 0.016e-2/2)
    T = 3600.0
    ρ = 0
    u_max(v) = 0.125
    u_min(v) = -0.25
    vᵢ = 2.0
    vf = 3.0

    prob = TrainProblem(track = track, resistance = myresistance, T = T, 
        umax = u_max, umin = u_min, ρ = ρ, vᵢ = vᵢ, vf = vf)
    solve!(prob)
end
