include("CircuitSim.jl")
using .CircuitSim

# Parameters of the circuit
C = 1e-3    # Capacitance (F)
R = 100.0     # Resistance (Ω)
L = 1e-3     # Inductance (H)
V₀ = 10.0   # Voltage amplitude (V)
f = 60.0     # Frequency (Hz)
ω = 2 * π * f # Angular frequency

# Creating the circuit 
circuit = Circuit(C, R, L, V₀, ω)
ε = define_voltage_cos(V₀, ω)

function circuit_ode(du, u, p, t)
    I_RL, I_C = u
    du[1] = (ε(t) - I_RL * circuit.R) / circuit.L       # dI_RL/dt (corrente através de R e L)
    du[2] = - circuit.C * ( ω^2 * ε(t))     # dI_C/dt (corrente através de C)
end

solve_ode_cos(circuit_ode, circuit, ε, 0.1)