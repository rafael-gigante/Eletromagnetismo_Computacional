import potential_initializer as PotentialInitializer
import potential_solver as PotentialSolver
import potential_plots as PotentialPlots

# Parameters
V1 = -1.0
V2 = 1.0
N = 31
Lx = 2.0
Ly = 2.0

# Initialize the potential box
V = PotentialInitializer.initialize_V_capacitor(V1, V2, N)

# Solve using Jacobi method
V_final = PotentialSolver.laplace_calculate(V)

# Plot results
PotentialPlots.equipotential_lines(V_final, Lx, Ly)
PotentialPlots.potential_surface(V_final, Lx, Ly)
PotentialPlots.electric_field(V_final, Lx, Ly)
