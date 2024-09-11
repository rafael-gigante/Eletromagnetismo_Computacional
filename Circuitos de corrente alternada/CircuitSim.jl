module CircuitSim

using DifferentialEquations, Plots

export  Circuit, define_voltage_sin, define_voltage_cos, solve_ode_sin, solve_ode_cos

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

function calculate_initial_conditions(circuit::Circuit, ε)
    I_c = circuit.C * circuit.ω * ε(pi/2)
    ϕ = atan((circuit.ω * circuit.L)/circuit.R)
    I_rl = - (circuit.V₀ * sin(ϕ))/sqrt(circuit.R^2 + (circuit.ω * circuit.L)^2)
    return u0 = [I_rl, I_c]
end

function solve_ode_sin(circuit_ode, circuit::Circuit, ε, tmax)
    u0 = calculate_initial_conditions(circuit, ε)
    tspan = (0.0, tmax)
    prob = ODEProblem(circuit_ode, u0, tspan)
    sol = solve(prob)
    
    time, I_c, I_rl = calculate_analytical_solutions_sin(circuit, ε, tmax)
    plot_results(sol,time, I_c, I_rl)
end

function solve_ode_cos(circuit_ode, circuit::Circuit, ε, tmax)
    u0 = calculate_initial_conditions(circuit, ε)
    tspan = (0.0, tmax)
    prob = ODEProblem(circuit_ode, u0, tspan)
    sol = solve(prob)
    
    time, I_c, I_rl = calculate_analytical_solutions_cos(circuit, ε, tmax)
    plot_results(sol,time, I_c, I_rl)
end

function I_c_analytical(circuit::Circuit, ε, t)
    return circuit.C * circuit.ω * ε(t + pi/2)
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
        push!(I_c, I_c_analytical(circuit, ε, i))
        push!(I_rl, I_rl_analytical_sin(circuit, i)) 
    end
    return time, I_c, I_rl
end

function calculate_analytical_solutions_cos(circuit::Circuit, ε, tmax)
    time = range(0, stop=tmax, length=200)

    I_c = []
    I_rl = []
    for i in time
        push!(I_c, I_c_analytical(circuit, ε, i))
        push!(I_rl, I_rl_analytical_cos(circuit, i)) 
    end
    return time, I_c, I_rl
end

function plot_results(sol, time, I_c, I_rl)
    pl = plot(sol, idxs=1,  xlabel = "Time (s)", ylabel = "Current (A)", title = "RL Branch", label="Solução Numérica",  color=:blue,  framestyle=:box, legend=:topright)
    scatter!(time, I_rl, label="Solução Analítica", color=:black, markersize=2)
    display(pl)

    pl = plot(sol, idxs=2, xlabel = "Time (s)", ylabel = "Current (A)", title = "Capacitor Branch", label="Solução Numérica", color=:blue,  framestyle=:box, legend=:topright)
    scatter!(time, I_c, label="Solução Analítica", color=:black, markersize=2)
    display(pl)
end

end