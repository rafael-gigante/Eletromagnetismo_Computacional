include("CircuitSim.jl")
using .CircuitSim

# Parameters of the circuit
C = 1e-6   # Capacitance (F)
R = 10.0     # Resistance (Ω)
L = 1e-3     # Inductance (H)
V₀ = 5.0   # Voltage amplitude (V)
f = 60.0     # Frequency (Hz)
ω = 2 * π * f # Angular frequency
ω₀ = sqrt((1/(L * C) - (R^2/L^2)))

# Creating the circuit 
circuit = Circuit(C, R, L, V₀, ω)
ε = define_voltage_sin(V₀, ω)

function circuit_ode(du, u, p, t)
    I_RL, I_C = u
    du[1] = (ε(t) - I_RL * circuit.R) / circuit.L       # dI_RL/dt (corrente através de R e L)
    du[2] = - circuit.C * ( ω^2 * ε(t))     # dI_C/dt (corrente através de C)
end

solve_ode_sin(circuit_ode, circuit, ε, 0.1)
println("Ressonance frequency:", ω₀) 
resonance_sin(circuit)