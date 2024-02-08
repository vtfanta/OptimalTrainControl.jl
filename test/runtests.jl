using OptimalTrainControl
using StaticArrays
using Test

@testset "Track utilities" begin
    t = Track(
        length = 1e3,
        altitude = 100.0,
        x_gradient = [0.0, 200.0, 500.0],
        gradient = [0.0, 1.0, 2.0],
        x_speedlimit = [0.0, 30.0, 600.0],
        speedlimit = [80., 85., 95.]
    )

    @test gradient(t, 94) ≈ 0
    @test gradient(t, 0) ≈ 0
    @test gradient(t, t.length) ≈ 2
    @test gradient(t, 250) ≈ 1
    @test gradient(t, 700) ≈ 2
    @test speedlimit(t, 20) ≈ 80
    @test speedlimit(t, 50) ≈ 85
    @test speedlimit(t, 700) ≈ 95
    @test speedlimit(t, 0) ≈ 80
    @test speedlimit(t, t.length) ≈ 95

    # check if track bounds checked correctly
    @test isvalidposition(t, -0.1) == false
    @test isvalidposition(t, t.length+1) == false
    @test isvalidposition(t, 10)

    @test altitude(t, 50) ≈ t.altitude
    @test altitude(t, 250) ≈ t.altitude + 50
    @test altitude(t, 600) ≈ t.altitude + 300 + 200
    @test altitude(t, t.length) ≈ t.altitude + 300 + 2*500

    @test all(segmentize!(t) .≈ [0,30,200,500,600])

    @test g(t, 250) ≈ -9.81*sqrt(2)/2
end

@testset "Train utilities" begin
    train = Train(
        v -> 1/v,
        v -> -1/v,
        (0.00675, 0.0, 0.00005),
        0.2
    )

    @test r(train, 15) ≈ 0.018
end

@testset "Time-optimal train control, flat track" begin
    train = Train(
        U̅ = v -> 3/v,
        U̲ = v -> -3/v,
        r = (6.75e-3, 0., 5e-5)
    )

    track = Track(
        length = 300.,
        altitude = 100.,
        x_gradient = [0.0],
        gradient = [2e-3]
    )

    timeprob = TOTCProblem(train, track)

    sol = solve(timeprob)
    T_end = sol.odesol.u[end][1]
    @test isapprox(T_end, 53.149; atol=0.1)
end

@testset "Time-optimal train control, wavy track" begin
    train = Train(
            U̅ = v -> 3/v,
            U̲ = v -> -3/v,
            r = (6.75e-3, 0., 5e-5)
        )

    track = Track(
        length = 300.,
        altitude = 100.,
        x_gradient = collect(0.:1.:300.),
        gradient = 5e-2*sin.(collect(0.:1.:300.)./20.)
    )

    timeprob = TOTCProblem(;train, track)

    sol = solve(timeprob)

    @test isapprox(sol.odesol[1,end], 54.30; atol = 0.1)
end

@testset "Energy-efficient train control, flat track" begin
    T = 800.
    
    train = Train(
            U̅ = v -> 3/v,
            U̲ = v -> -3/v,
            r = (6.75e-3, 0., 5e-5)
        )

    track = Track(length = 5e3)

    prob = EETCProblem(; T, train, track)
    sol = solve(prob)
    
    @test isapprox(sol.odesol[1,end], T; atol = 5.)
    @test isapprox(sol.odesol[2,end], 1.; atol = 0.1)
end