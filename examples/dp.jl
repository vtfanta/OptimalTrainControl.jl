using BasicInterpolators
using Plots
using OptimalTrainControl

# From HW5 of ORR

function twostateDP(dt, prob::TrainProblem)
    function f(x, u, t)
        # Euler's method discretization
        if x[2] < 0
            @show x[2]
        end
        xnext = zeros(2)
        xnext[1] = x[1] + dt * 1/x[2]
        xnext[2] = x[2] + dt * (u - resistance(prob.resistance, x[2]) + 
            getgradientacceleration(prob.track, t)) / x[2]
        return xnext
    end

    # running cost
    L(x, u) = abs(u)
    # terminal cost
    ϕ(x) = 1000000 * (abs2(x[2] - prob.vf) + abs2(x[1] - prob.T))

    # amount of discretization steps in states
    stateK = 15
    xb = zeros(2, stateK)
    xb[1,:] = collect(range(0, prob.T+10; length = stateK))
    xb[2,:] = collect(range(0.1, 25; length = stateK))

    # amount of discretization steps in control
    controlK = 35
    ub = range(prob.umin(0), prob.umax(0); length = controlK)

    # Time discretization
    t = range(0, finish(prob.track); step = dt)

    uopt = zeros(length(t) - 1, stateK, stateK)
    J = zeros(length(t), stateK, stateK)

    # J in final step calculation
    for i = 1:stateK, j = 1:stateK
        J[end, i, j] = ϕ([xb[1,i], xb[2,j]])
    end

    for k ∈ length(t) - 1:-1:1
        @show k
        for x1idx ∈ 1:stateK, x2idx ∈ 1:stateK
            Jk_min = Inf
            uk_opt = nothing
            for uidx ∈ 1:controlK
                if !(prob.umin(xb[2,x2idx]) ≤ ub[uidx] ≤ prob.umax(xb[2,x2idx]))
                    continue
                end

                x_next = f([xb[1, x1idx], xb[2, x2idx]], ub[uidx], t[k])
                if !(xb[1,1] ≤ x_next[1] ≤ xb[1,end]) || !(xb[2,1] ≤ x_next[2] ≤ xb[2,end])
                    continue
                end

                itp = BilinearInterpolator(xb[1,:], xb[2,:], J[k+1,:,:], NoBoundaries()) 
                
                Jk_next = itp(x_next...)

                Jk = L([xb[1, x1idx], xb[2, x2idx]], ub[uidx]) + Jk_next
                if Jk < Jk_min
                    Jk_min = Jk
                    uk_opt = ub[uidx]
                end
            end

            if !isinf(Jk_min)
                J[k,x1idx,x2idx] = Jk_min   
                uopt[k, x1idx, x2idx] = uk_opt
            else
                J[k,x1idx,x2idx] = Inf
            end
        end
    end

    # simulation
    x = zeros(2, length(t))
    u = zeros(length(t) - 1)

    x[:,1] = [0, prob.vᵢ]

    for k ∈ 2:length(t)
        itp = BilinearInterpolator(xb[1,:], xb[2,:], uopt[k-1,:,:], NoBoundaries())
        u[k-1] = itp(x[:,k-1]...)
        x[:,k] = f(x[:,k-1],u[k-1],t[k])
    end

    return t, x, u
end

function hw5_fantavit()
    function f(x, u)
        x - 0.4x^2 + u
    end

    L(x, u) = abs(x)
    ϕ(x) = abs2(x - 0.5)

    xb = -1:0.1:1
    ub = -0.5:0.1:0.5

    t = 0:5

    uopt = zeros(length(t) - 1, length(xb))
    J = zeros(length(t), length(xb))

    J[end, :] = ϕ.(xb)

    for k ∈ length(t) - 1:-1:1
        for i ∈ eachindex(xb)
            Jk_min = Inf
            uk_opt = nothing
            for j ∈ eachindex(ub)
                x_next = f(xb[i], ub[j])
                if x_next < xb[1] || x_next > xb[end]
                    continue
                end

                itp = LinearInterpolator(xb, J[k+1,:], NoBoundaries())
                Jk_next = itp(x_next)

                Jk = L(xb[i], ub[j]) + Jk_next

                if Jk < Jk_min
                    Jk_min = Jk
                    uk_opt = ub[j]
                end
            end
            J[k,i] = Jk_min
            uopt[k,i] = uk_opt
        end
    end

    x = zeros(length(t))
    x[1] = 0.1
    u = zeros(length(t)-1)

    for i = 2:length(t)
        itp = LinearInterpolator(xb, uopt[i-1,:])
        u[i-1] = itp(x[i-1])
        x[i] = f(x[i-1], u[i-1])
    end

    t, x, u
end

# t, x, u = hw5_fantavit()
track = FlatTrack(5e3)
prob = TrainProblem(;track, T = 1000, umax = v -> 0.125, umin = v -> -0.25)
t, x, u = twostateDP(20, prob)
l = @layout [a;b;c]
px1 = plot(t, x[1,:], label = "x1")
px2 = plot(t, x[2,:], label = "x2")
pu = plot(t[1:end-1], u, label = "u")
plot(px1, px2, pu, layout = l)
