import potential_initializer as PotentialInitializer
import potential_solver as PotentialSolver
import potential_plots as PotentialPlots

# Parameters
V1 = -1.0
V2 = 1.0
N = 31
Lx = 2.0
Ly = 2.0
methods = ["Jacobi", "GS", "SOR"]
method = methods[2]

# Initialize the potential box
V, ρ = PotentialInitializer.initialize_V_box_charge(N)

dim = len(V.shape)

if dim == 2:
    # Solve using Jacobi method
    V_final = PotentialSolver.laplace_calculate(V, ρ, method)

    # Plot results
    PotentialPlots.equipotential_lines(V_final, Lx, Ly)
    PotentialPlots.potential_surface(V_final, Lx, Ly)
    PotentialPlots.electric_field(V_final, Lx, Ly)
elif dim ==3:
    for _ in range(100):
        V_new, ΔV = PotentialSolver.update_V(V, ρ, method)
        V, ΔV = PotentialSolver.update_V(V_new, ρ, method)
    V_final = V

    z0 = int((N-1)/2) + 1
    V_xy = V_final[:, :, z0]

    # Plot results
    PotentialPlots.equipotential_lines(V_xy, Lx, Ly)
    PotentialPlots.potential_surface(V_xy, Lx, Ly)
    PotentialPlots.electric_field(V_xy, Lx, Ly)


        

    
