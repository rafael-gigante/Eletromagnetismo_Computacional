module CircuitSim

using DifferentialEquations, Plots

export  Circuit, define_voltage_sin, define_voltage_cos, solve_ode_sin, solve_ode_cos, resonance_sin

struct Circuit
    C::Float64 
    R::Float64 
    L::Float64 
    V₀::Float64 
    ω::Float64
end

function define_voltage_sin(V₀, ω)
    return ε(t) = V₀ * sin(ω * t)
end

function define_voltage_cos(V₀, ω)
    return ε(t) = V₀ * cos(ω * t)
end

function calculate_initial_conditions_sin(circuit::Circuit, ε)
    I_c = circuit.C * circuit.ω * circuit.V₀
    ϕ = atan((circuit.ω * circuit.L)/circuit.R)
    I_rl = - (circuit.V₀ * sin(ϕ))/sqrt(circuit.R^2 + (circuit.ω * circuit.L)^2)
    return u0 = [I_rl, I_c]
end

function calculate_initial_conditions_cos(circuit::Circuit, ε)
    I_c = 0
    ϕ = atan((circuit.ω * circuit.L)/circuit.R)
    I_rl = (circuit.V₀ * cos(ϕ))/sqrt(circuit.R^2 + (circuit.ω * circuit.L)^2)
    return u0 = [I_rl, I_c]
end

function solve_ode_sin(circuit_ode, circuit::Circuit, ε, tmax)
    u0 = calculate_initial_conditions_sin(circuit, ε)
    tspan = (0.0, tmax)
    prob = ODEProblem(circuit_ode, u0, tspan)
    sol = solve(prob)
    
    time, I_c, I_rl = calculate_analytical_solutions_sin(circuit, ε, tmax)
    plot_results(circuit, sol,time, I_c, I_rl, ε)
end

function solve_ode_cos(circuit_ode, circuit::Circuit, ε, tmax)
    u0 = calculate_initial_conditions_cos(circuit, ε)
    tspan = (0.0, tmax)
    prob = ODEProblem(circuit_ode, u0, tspan)
    sol = solve(prob)
    
    time, I_c, I_rl = calculate_analytical_solutions_cos(circuit, ε, tmax)
    plot_results(circuit, sol, time, I_c, I_rl, ε)
end

function I_c_analytical_sin(circuit::Circuit, ε, t)
    return circuit.C * circuit.ω * circuit.V₀ * cos(circuit.ω * t)
end

function I_c_analytical_cos(circuit::Circuit, ε, t)
    return - circuit.C * circuit.ω * circuit.V₀ * sin(circuit.ω * t)
end

function I_rl_analytical_sin(circuit::Circuit, t)
    ϕ = atan((circuit.ω * circuit.L)/circuit.R)
    return (circuit.V₀ * sin(circuit.ω * t - ϕ))/sqrt(circuit.R^2 + (circuit.ω * circuit.L)^2)
end

function I_rl_analytical_cos(circuit::Circuit, t)
    ϕ = atan((circuit.ω * circuit.L)/circuit.R)
    return (circuit.V₀ * cos(circuit.ω * t - ϕ))/sqrt(circuit.R^2 + (circuit.ω * circuit.L)^2)
end

function calculate_analytical_solutions_sin(circuit::Circuit, ε, tmax)
    time = range(0, stop=tmax, length=200)

    I_c = []
    I_rl = []
    for i in time
        push!(I_c, I_c_analytical_sin(circuit, ε, i))
        push!(I_rl, I_rl_analytical_sin(circuit, i)) 
    end
    return time, I_c, I_rl
end

function calculate_analytical_solutions_cos(circuit::Circuit, ε, tmax)
    time = range(0, stop=tmax, length=200)

    I_c = []
    I_rl = []
    for i in time
        push!(I_c, I_c_analytical_cos(circuit, ε, i))
        push!(I_rl, I_rl_analytical_cos(circuit, i)) 
    end
    return time, I_c, I_rl
end

function plot_results(circuit, sol, time, I_c, I_rl, ε)
    V = []
    for i in time
        push!(V, ε(i)/circuit.V₀)
    end

    t_sol = sol.t
    I1 = [u[1] for u in sol.u]  
    I2 = [u[2] for u in sol.u]  

    I_max_rl = maximum(I_rl)
    I_max_c = maximum(I_c)
    τ = circuit.L / circuit.R

    I_tot = I1 .+ I2
    I_max_tot = maximum(I_tot)

    pl = plot(xlabel = "t/τ", ylabel = "Relative Amplitude", title = "RL Branch with R = $(circuit.R) and L = $(circuit.L)", framestyle=:box, legend=:topright, minorgrid=true)
    plot!(t_sol./τ, I1./I_max_rl, label="Numerical Solution",  color=:blue, linewidth=2)
    scatter!(time./τ, I_rl./I_max_rl, label="Analytical Solution", color=:black, markersize=3)
    plot!(time./τ, V, label="AC voltage",  color=:red, linewidth=1)
    display(pl)
    savefig(pl, "Circuitos de corrente alternada/figures/rl_current.png")

    pl = plot(xlabel = "t/τ", ylabel = "Relative Amplitude", title = "Capacitor Branch with C = $(circuit.C)",  framestyle=:box, legend=:topright, minorgrid=true)
    plot!(t_sol./τ, I2./I_max_c, label="Numerical Solution", color=:blue, linewidth=2)
    scatter!(time./τ, I_c./I_max_c, label="Analytical Solution", color=:black, markersize=3)
    plot!(time./τ, V, label="AC voltage",  color=:red, linewidth=1)
    display(pl)
    savefig(pl, "Circuitos de corrente alternada/figures/c_current.png")

    pl = plot(xlabel = "t/τ", ylabel = "I/Iₘₐₓ", title = "RL-C circuit (Numerical results)",  framestyle=:box, legend=:topright, minorgrid=true)
    plot!(t_sol./τ, I_tot./I_max_tot, label="Total current", color=:red, linewidth=2)
    plot!(t_sol./τ, I1./I_max_tot, label="RL current",  color=:green, linewidth=1)
    plot!(t_sol./τ, I2./I_max_tot, label="C current", color=:blue, linewidth=1)
    display(pl)
    savefig(pl, "Circuitos de corrente alternada/figures/total_current.png")
end

function resonance_sin(circuit::Circuit)
    ω₀ = sqrt((1/(circuit.L * circuit.C) - (circuit.R^2/circuit.L^2)))
    ω_range = range(0, 2*ω₀, length=6000)
    I_amplitude = []
    for w in ω_range
        ε = define_voltage_sin(circuit.V₀, w)

        time, I_c, I_rl = calculate_analytical_solutions_sin(circuit, ε, 0.1)
        I_tot = sqrt.(I_rl.^2 .+ I_c.^2)
        I_max = maximum(I_tot)
        push!(I_amplitude, I_max)
    end  
    pl=plot(ω_range, I_amplitude, xlabel="Angular Frequency (rad/s)", ylabel="Current Amplitude (A)",
     label="Current Amplitude", lw=2, legend=:topright)
    vline!([ω₀], label="Resonance Frequency", lw=2, color=:red)
    title!("Current Amplitude vs Angular Frequency")
    ylims!(0.0, 1.0)
    display(pl)
end

function resonance_ode_sin(circuit_ode, circuit::Circuit, ε, tmax)
    u0 = calculate_initial_conditions_sin(circuit, ε)
    tspan = (0.0, tmax)
    prob = ODEProblem(circuit_ode, u0, tspan)
    sol = solve(prob)
    I1 = [u[1] for u in sol.u]  
    I2 = [u[2] for u in sol.u]  
    I_max_rl = maximum(I1)
    I_max_c = maximum(I2)
    I_max = I_max_rl + I_max_c
    return I_max
end

end