import numpy as np

def update_V(V, ρ, method):
    """
    Function that receives a potential V and a charge distribution and aplies chosen method to it
    and calculates the absolute difference between the new and old potentials.

    Returns:
    - `V_new`: Array with the new potential calculate by the Jacobi method,
    - 'ΔV': Float of the absolute difference between the potential received and the calculated.
    """
    if (method == "Jacobi"):
        V_new, ΔV = jacobi_method(V, ρ)
    elif (method == "GS"):
        V_new, ΔV = gauss_seidel_method(V, ρ)
    elif (method == "SOR"):
        V_new, ΔV = sor_method(V, ρ)
        
    return V_new, ΔV

def jacobi_method(V, ρ):

    # Total (cumulative) change in V during this update loop
    ΔV = 0

    # Dimension of the potential
    dim = len(V.shape)

    # Number of sites along a direction
    N = V.shape[0]

    # Step size
    dx = 2.0 / (N - 1)

    V_new = np.copy(V)
    # Loop through all values (i,j, k) except the boundary, where the values are fixed
    if dim == 2:
        for i in range(1, N-1):
            for j in range(1, N-1):
                if (abs(V[i, j]) != 1):
                    V_new[i, j] = (V[i+1, j] + V[i-1, j] +
                                V[i, j+1] + V[i, j-1] + ρ[i, j]*(dx**2)) / 4
                    ΔV += abs(V_new[i, j] - V[i, j])
    elif dim == 3:
        for i in range(1, N-1):
            for j in range(1, N-1):
                for k in range(1, N-1):
                    V_new[i, j, k] = (V[i+1, j, k] + V[i-1, j, k] +
                                V[i, j+1, k] + V[i, j-1, k] +
                                V[i, j, k+1] + V[i, j, k-1] + ρ[i, j, k]*(dx**2)) / 6
                    ΔV += abs(V_new[i, j, k] - V[i, j, k])
    return V_new, ΔV

def gauss_seidel_method(V, ρ):

    # Total (cumulative) change in V during this update loop
    ΔV = 0

    # Dimension of the potential
    dim = len(V.shape)

    # Number of sites along a direction
    N = V.shape[0]

    # Step size
    dx = 2.0 / (N - 1)

    V_new = np.copy(V)
    # Loop through all values (i,j, k) except the boundary, where the values are fixed
    if dim == 2:
        for i in range(1, N-1):
            for j in range(1, N-1):
                if (abs(V[i, j]) != 1):
                    V_new[i, j] = (V[i+1, j] + V_new[i-1, j] +
                                V[i, j+1] + V_new[i, j-1] + ρ[i, j]*(dx**2)) / 4
                    ΔV += abs(V_new[i, j] - V[i, j])
    elif dim == 3:
        for i in range(1, N-1):
            for j in range(1, N-1):
                for k in range(1, N-1):
                    V_new[i, j, k] = (V[i+1, j, k] + V_new[i-1, j, k] +
                                V[i, j+1, k] + V_new[i, j-1, k] +
                                V[i, j, k+1] + V_new[i, j, k-1] + ρ[i, j, k]*(dx**2)) / 6
                    ΔV += abs(V_new[i, j, k] - V[i, j, k])
    return V_new, ΔV

def sor_method(V, ρ):

    # Number of sites along a direction
    N = V.shape[0]

    # Dimension of the potential
    dim = len(V.shape)

    α = 1.5

    V_new = np.copy(V)
    V_new, ΔV = gauss_seidel_method(V, ρ)

    V_new = α * (V_new - V) + V

    # Total (cumulative) change in V during this update loop
    ΔV = 0

    # Loop through all values (i,j, k) except the boundary, where the values are fixed
    if dim == 2:
        for i in range(1, N-1):
            for j in range(1, N-1):
                ΔV += abs(V_new[i, j] - V[i, j])
    elif dim == 3:
        for i in range(1, N-1):
            for j in range(1, N-1):
                for k in range(1, N-1):
                    ΔV += abs(V_new[i, j, k] - V[i, j, k])
    return V_new, ΔV


def laplace_calculate(V, ρ, method):
    """
    Function that aplies the relaxation method to a initial potential V
    until the new potential is within the tolerance range.

    Returns:
    - `V`: Array with final potential that is within the tolerance range
    """

    # Define the error tolerance
    N = V.shape[0]
    ε = 10**(-6) * N**2

    # Perform 10 iterations initially
    for _ in range(10):
        V_new, ΔV = update_V(V, ρ, method)
        V, ΔV = update_V(V_new, ρ, method)

    # Continue iterating until ΔV is less than ε
    while ΔV < ε:
        V_new, ΔV = update_V(V, ρ, method)
        V, ΔV = update_V(V_new, ρ, method)

    return V
